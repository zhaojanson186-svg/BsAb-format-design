import streamlit as st
import re
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ==========================================
# 核心底层类：BsAb Assembler
# ==========================================
class BsAbAssembler:
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

    def __init__(self, fc_silence, linker, lc_type):
        self.fc_silence = fc_silence
        self.linker = linker
        self.lc_type = lc_type

    def clean_sequence(self, seq):
        if not seq: return ""
        # 移除非字母字符并转大写，同时移除可能导致 Biopython 报错的非常规氨基酸占位符 (如 X, B, Z)
        clean_seq = re.sub(r'[^a-zA-Z]', '', seq).upper()
        return re.sub(r'[XBZ]', '', clean_seq)

    def analyze_protein(self, seq):
        try:
            analyzer = ProteinAnalysis(seq)
            mw = analyzer.molecular_weight()
            pi = analyzer.isoelectric_point()
            return round(mw / 1000, 2), round(pi, 2)
        except Exception:
            return 0.0, 0.0

    def assemble(self, taa_vh, taa_vl, cd3_vh, cd3_vl, project_name):
        taa_vh = self.clean_sequence(taa_vh)
        taa_vl = self.clean_sequence(taa_vl)
        cd3_vh = self.clean_sequence(cd3_vh)
        cd3_vl = self.clean_sequence(cd3_vl)

        link = self.TEMPLATES["linker"][self.linker]
        ch1_hinge = self.TEMPLATES["ch1_hinge_front"]["IgG1_WT"]
        fc_knob = self.TEMPLATES["fc_knob"][self.fc_silence]
        fc_hole = self.TEMPLATES["fc_hole_r"][self.fc_silence]
        ck = self.TEMPLATES["light_constant"][self.lc_type]

        chain_h1 = taa_vh + ch1_hinge + link + cd3_vl + link + cd3_vh + fc_knob
        chain_h2 = taa_vh + fc_hole
        chain_l1 = taa_vl + ck

        chains = {
            f"{project_name}_H1 (scFv-Fc_Knob)": chain_h1,
            f"{project_name}_H2 (Fab-Fc_Hole)": chain_h2,
            f"{project_name}_L1 (Common_LC)": chain_l1
        }
        
        results = []
        for name, seq in chains.items():
            mw_kda, pi = self.analyze_protein(seq)
            results.append({
                "Chain": name,
                "Sequence": seq,
                "Length (AA)": len(seq),
                "MW (kDa)": mw_kda,
                "pI": pi
            })
        return results

# ==========================================
# Streamlit 前端 UI 布局
# ==========================================
st.set_page_config(page_title="Next-Gen TCE Assembler", page_icon="🧬", layout="wide")

st.title("🧬 新一代 TCE 双抗 (2+1 构型) 组装与评估平台")
st.markdown("只需输入单抗可变区序列，系统将自动进行**数据清洗**、**架构拼装**及**理化性质计算**。底层已内置 LALA-GA 沉默与 KiH (H435R) 纯化优化设计。")

# --- 侧边栏：模板与参数配置 ---
with st.sidebar:
    st.header("⚙️ 模板配置 (Templates)")
    st.markdown("---")
    project_name = st.text_input("📦 项目/克隆名称", value="BMK_B10")
    fc_option = st.selectbox("🛡️ Fc 沉默突变", ["LALA_GA", "LALA_PG"])
    linker_option = st.selectbox("🔗 scFv 柔性接头", ["G4S_3", "Rigid_18"])
    lc_option = st.selectbox("⛓️ 公共轻链类型", ["Kappa", "Lambda"])
    
    st.markdown("---")
    st.info("**CMC 提示:**\n- H1/H2 链的 pI 差异可指导阳离子交换层析 (CEX) 策略。\n- 理论 MW 将作为 LC-MS 质谱鉴定的金标准。")

# --- 主界面：序列输入区 ---
col1, col2 = st.columns(2)

with col1:
    st.subheader("🎯 TAA (肿瘤靶向臂)")
    taa_vh_input = st.text_area("TAA VH 序列", height=150, placeholder="粘贴 VH 序列 (包含或不包含换行、空格均可)...")
    taa_vl_input = st.text_area("TAA VL 序列", height=150, placeholder="粘贴 VL 序列 (包含或不包含换行、空格均可)...")

with col2:
    st.subheader("⚔️ CD3 (T细胞接合臂)")
    cd3_vh_input = st.text_area("CD3 VH 序列", height=150, placeholder="粘贴 VH 序列 (包含或不包含换行、空格均可)...")
    cd3_vl_input = st.text_area("CD3 VL 序列", height=150, placeholder="粘贴 VL 序列 (包含或不包含换行、空格均可)...")

# --- 运行逻辑 ---
if st.button("🚀 开始组装双抗序列", type="primary", use_container_width=True):
    if not all([taa_vh_input, taa_vl_input, cd3_vh_input, cd3_vl_input]):
        st.warning("⚠️ 请填满所有四个可变区 (VH/VL) 序列。")
    else:
        with st.spinner('正在进行底层拼装与理化计算...'):
            assembler = BsAbAssembler(fc_silence=fc_option, linker=linker_option, lc_type=lc_option)
            results = assembler.assemble(taa_vh_input, taa_vl_input, cd3_vh_input, cd3_vl_input, project_name)
            
            st.success("✅ 组装成功！所有非法字符及换行已自动清洗。")
            
            # 1. 渲染理化性质表格
            st.subheader("📊 理化性质分析 (Physicochemical Parameters)")
            df = pd.DataFrame(results)[["Chain", "Length (AA)", "MW (kDa)", "pI"]]
            st.dataframe(df, use_container_width=True, hide_index=True)
            
            # 2. 渲染 FASTA 格式结果
            st.subheader("📝 组装结果 (FASTA)")
            fasta_output = ""
            for item in results:
                fasta_output += f">{item['Chain']} | Length: {item['Length (AA)']} | MW: {item['MW (kDa)']}kDa | pI: {item['pI']}\n"
                seq = item['Sequence']
                # 每 80 个氨基酸换行，符合 FASTA 规范
                for i in range(0, len(seq), 80):
                    fasta_output += seq[i:i+80] + "\n"
                fasta_output += "\n"
            
            st.code(fasta_output, language="fasta")
