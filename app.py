import streamlit as st
import pandas as pd
import re
from io import BytesIO
from docx import Document
from docx.shared import RGBColor, Pt
from docx.enum.text import WD_COLOR_INDEX
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ==========================================
# 1. 底层序列元件与模板库
# ==========================================
TEMPLATES = {
    "linker": {
        "G4S_3": "GGGGSGGGGSGGGGS",
        "Rigid_18": "KPGSGKPGSGKPGSGKPG"
    },
    "ch1_hinge_front": {
        "IgG1_WT": "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC"
    },
    "fc_knob": {
        "LALA_GA": "EPKSSDKTHTCPPCPAPEAAGAPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPCRDELTKNQVSLWCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK",
        "LALA_PG": "EPKSSDKTHTCPPCPAPEAAGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALGAPIEKTISKAKGQPREPQVYTLPPCRDELTKNQVSLWCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    },
    "fc_hole_r": {
        "LALA_GA": "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPEAAGAPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVCTLPPSRDELTKNQVSLSCAVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLVSKLTVDKSRWQQGNVFSCSVMHEALHNRYTQKSLSLSPGK",
        "LALA_PG": "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPEAAGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALGAPIEKTISKAKGQPREPQVCTLPPSRDELTKNQVSLSCAVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLVSKLTVDKSRWQQGNVFSCSVMHEALHNRYTQKSLSLSPGK"
    },
    "light_constant": {
        "Kappa": "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
        "Lambda": "GQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    }
}

# ==========================================
# 2. 颜色配置字典 (用于 Word 导出)
# ==========================================
COLOR_MAP = {
    'TAA_VH': RGBColor(31, 78, 121),    # 深蓝
    'TAA_VL': RGBColor(0, 176, 240),    # 浅蓝
    'CD3_VH': RGBColor(192, 0, 0),      # 深红
    'CD3_VL': RGBColor(255, 0, 0),      # 亮红
    'Linker': RGBColor(56, 87, 35),     # 森林绿
    'Constant': RGBColor(127, 127, 127),# 灰色 (CH1, Ck)
    'Fc': RGBColor(197, 90, 17)         # 橙色 (Knob/Hole)
}

# ==========================================
# 3. 核心功能函数
# ==========================================
def parse_fasta(text):
    """解析批量 FASTA 输入，通过名字分组提取 VH 和 VL"""
    clones = {}
    current_name = None
    current_chain = None
    
    for line in text.split('\n'):
        line = line.strip()
        if not line: continue
        
        if line.startswith('>'):
            # 期望格式: >CloneName_VH 或 >CloneName_VL
            header = line[1:].strip()
            parts = header.split('_')
            if len(parts) >= 2:
                current_chain = parts[-1].upper() # 提取最后一部分 (VH/VL)
                current_name = '_'.join(parts[:-1]) # 剩余部分作为克隆名
                if current_name not in clones:
                    clones[current_name] = {}
            else:
                current_name = header
                current_chain = 'UNKNOWN'
        else:
            if current_name and current_chain in ['VH', 'VL']:
                # 清洗序列
                seq = re.sub(r'[^a-zA-Z]', '', line).upper()
                if current_chain not in clones[current_name]:
                    clones[current_name][current_chain] = ""
                clones[current_name][current_chain] += seq
                
    # 验证完整性 (只保留同时拥有 VH 和 VL 的克隆)
    valid_clones = {k: v for k, v in clones.items() if 'VH' in v and 'VL' in v}
    return valid_clones

def assemble_segments(taa_vh, taa_vl, cd3_vh, cd3_vl, format_type, fc_silence, linker, lc_type):
    """根据选择的构型，返回分段的序列对象列表，以便后续上色"""
    chains = {}
    
    # 调取模板
    link_seq = TEMPLATES["linker"][linker]
    ch1_hinge = TEMPLATES["ch1_hinge_front"]["IgG1_WT"]
    fc_knob = TEMPLATES["fc_knob"][fc_silence]
    fc_hole = TEMPLATES["fc_hole_r"][fc_silence]
    ck = TEMPLATES["light_constant"][lc_type]
    
    if format_type == "2+1 构型 (TAA双臂 + CD3单臂)":
        chains['H1 (scFv-Fc_Knob)'] = [
            ('TAA_VH', taa_vh), ('Constant', ch1_hinge), ('Linker', link_seq),
            ('CD3_VL', cd3_vl), ('Linker', link_seq), ('CD3_VH', cd3_vh), ('Fc', fc_knob)
        ]
        chains['H2 (Fab-Fc_Hole)'] = [
            ('TAA_VH', taa_vh), ('Fc', fc_hole)  # fc_hole_r 模板内已经包含了 CH1 和 Hinge
        ]
        chains['L1 (Common_LC)'] = [
            ('TAA_VL', taa_vl), ('Constant', ck)
        ]
        
    elif format_type == "1+1 构型 (TAA单臂 + CD3单臂)":
        chains['H1 (CD3_scFv-Fc_Knob)'] = [
            ('CD3_VL', cd3_vl), ('Linker', link_seq), ('CD3_VH', cd3_vh), ('Fc', fc_knob)
        ]
        chains['H2 (TAA_Fab-Fc_Hole)'] = [
            ('TAA_VH', taa_vh), ('Fc', fc_hole)
        ]
        chains['L1 (TAA_LC)'] = [
            ('TAA_VL', taa_vl), ('Constant', ck)
        ]
        
    return chains

def generate_word_document(combinations_data):
    """生成带颜色标记的 Word 文档"""
    doc = Document()
    doc.styles['Normal'].font.name = 'Courier New'
    doc.styles['Normal'].font.size = Pt(10)
    
    title = doc.add_heading('HTS Bispecific Antibody Sequence Report', 0)
    
    # 图例
    legend = doc.add_paragraph()
    legend.add_run("Color Legend: ").bold = True
    for key, color in COLOR_MAP.items():
        run = legend.add_run(f" [{key}] ")
        run.font.color.rgb = color
        run.bold = True
        
    for combo in combinations_data:
        doc.add_heading(combo['Name'], level=1)
        doc.add_paragraph(f"Combinatorial Pair: TAA [{combo['TAA_Clone']}] x CD3 [{combo['CD3_Clone']}]")
        
        for chain_name, segments in combo['Chains'].items():
            p = doc.add_paragraph()
            p.add_run(f">{combo['Name']}_{chain_name}\n").bold = True
            
            # 拼装颜色块
            for seg_type, seq in segments:
                run = p.add_run(seq)
                if seg_type in COLOR_MAP:
                    run.font.color.rgb = COLOR_MAP[seg_type]
                    
    # 保存到内存
    bio = BytesIO()
    doc.save(bio)
    return bio.getvalue()

# ==========================================
# 4. Streamlit UI 布局
# ==========================================
st.set_page_config(page_title="TCE HTS Combinatorial Platform", page_icon="🔬", layout="wide")

st.title("🔬 工业级 TCE 双抗高通量组合与生成平台")
st.markdown("通过批量输入候选克隆 FASTA 序列，自动执行**交叉组合（Cross-pairing）**、自增命名，并生成含**模块可视化着色**的研发级 Word 报告。")

# --- 侧边栏 ---
with st.sidebar:
    st.header("⚙️ 平台构建参数")
    st.markdown("---")
    base_name = st.text_input("📦 基础项目命名 (Base Name)", value="BMK")
    start_idx = st.number_input("🔢 起始编号", min_value=1, value=1)
    
    st.markdown("---")
    format_option = st.radio("📐 双抗构型 (Format)", ["2+1 构型 (TAA双臂 + CD3单臂)", "1+1 构型 (TAA单臂 + CD3单臂)"])
    fc_option = st.selectbox("🛡️ Fc 沉默底盘", ["LALA_GA", "LALA_PG"])
    linker_option = st.selectbox("🔗 scFv 接头类型", ["G4S_3", "Rigid_18"])
    lc_option = st.selectbox("⛓️ 恒定区轻链类型", ["Kappa", "Lambda"])

# --- 主界面：批量输入区 ---
st.info("💡 **格式提示**: 请确保输入的 FASTA 包含克隆名及对应链的后缀标识。例如：`>TAA01_VH` 和 `>TAA01_VL`。系统将自动按名字进行组装匹配。")

col1, col2 = st.columns(2)
with col1:
    st.subheader("🎯 TAA 克隆文库 (Batch Input)")
    taa_fasta = st.text_area("输入 TAA 序列", height=250, placeholder=">CloneA_VH\nEVQL...\n>CloneA_VL\nDIQM...")
with col2:
    st.subheader("⚔️ CD3 克隆文库 (Batch Input)")
    cd3_fasta = st.text_area("输入 CD3 序列", height=250, placeholder=">SP34_VH\nEVQL...\n>SP34_VL\nQTVV...")

# --- 核心运行逻辑 ---
if st.button("🚀 批量交叉组合并生成文件", type="primary", use_container_width=True):
    if not taa_fasta or not cd3_fasta:
        st.warning("⚠️ 必须同时提供 TAA 和 CD3 的批量序列！")
    else:
        with st.spinner("正在解析序列并执行组合拼装..."):
            taa_clones = parse_fasta(taa_fasta)
            cd3_clones = parse_fasta(cd3_fasta)
            
            if not taa_clones:
                st.error("❌ 无法解析 TAA 序列，请检查是否按 `>克隆名_VH` 和 `>克隆名_VL` 格式输入。")
            elif not cd3_clones:
                st.error("❌ 无法解析 CD3 序列，请检查格式。")
            else:
                # 执行交叉组合
                combinations = []
                counter = start_idx
                
                for taa_name, taa_seqs in taa_clones.items():
                    for cd3_name, cd3_seqs in cd3_clones.items():
                        combo_name = f"{base_name}-{counter:03d}"
                        counter += 1
                        
                        # 获取分段序列对象
                        chains = assemble_segments(
                            taa_vh=taa_seqs['VH'], taa_vl=taa_seqs['VL'],
                            cd3_vh=cd3_seqs['VH'], cd3_vl=cd3_seqs['VL'],
                            format_type=format_option, fc_silence=fc_option, 
                            linker=linker_option, lc_type=lc_option
                        )
                        
                        combinations.append({
                            'Name': combo_name,
                            'TAA_Clone': taa_name,
                            'CD3_Clone': cd3_name,
                            'Chains': chains
                        })
                
                st.success(f"✅ 成功生成 {len(combinations)} 种双抗组合！ (TAA克隆数: {len(taa_clones)} × CD3克隆数: {len(cd3_clones)})")
                
                # 展现总览表格
                overview_data = [{"分子编号": c['Name'], "TAA 来源": c['TAA_Clone'], "CD3 来源": c['CD3_Clone']} for c in combinations]
                st.dataframe(pd.DataFrame(overview_data), use_container_width=True)
                
                # 生成 Word 文档并提供下载
                word_file = generate_word_document(combinations)
                file_name = f"{base_name}_HTS_Library_{format_option.split(' ')[0]}.docx"
                
                st.download_button(
                    label="💾 下载着色序列报告 (Word 文档)",
                    data=word_file,
                    file_name=file_name,
                    mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                    type="primary"
                )
