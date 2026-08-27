import streamlit as st
import requests
import urllib.parse
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, rdMolDescriptors

# --- 1. PAGE CONFIGURATION ---
st.set_page_config(
    page_title="FormuAI Computational Engine",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- 2. CUSTOM DARK METRIC & DASHBOARD STYLING ---
st.markdown("""
<style>
    .stApp {
        background-color: #0E1117;
        color: #E0E6ED;
    }
    .main-header {
        background: linear-gradient(135deg, #1E293B 0%, #0F172A 100%);
        border: 1px solid #334155;
        border-radius: 12px;
        padding: 20px;
        margin-bottom: 20px;
    }
    .main-header h1 {
        color: #38BDF8;
        font-size: 1.8rem;
        margin-bottom: 4px;
    }
    .metric-card {
        background: #1E293B;
        border: 1px solid #334155;
        border-radius: 8px;
        padding: 16px;
        text-align: center;
    }
    .metric-val {
        font-size: 1.6rem;
        font-weight: 700;
        color: #38BDF8;
    }
    .metric-lbl {
        font-size: 0.8rem;
        color: #94A3B8;
        text-transform: uppercase;
    }
    .footer-container {
        margin-top: 40px;
        padding: 16px;
        border-top: 1px solid #334155;
        background-color: #0F172A;
        border-radius: 8px;
        color: #94A3B8;
        font-size: 0.85rem;
    }
</style>
""", unsafe_allow_html=True)

# --- 3. SESSION STATE INITIALIZATION ---
if "compound_name" not in st.session_state:
    st.session_state.compound_name = "diazepam"
if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3"

# --- 4. PUBCHEM API ENGINE WITH TIMEOUT & HEADERS ---
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_pubchem_data(query_name):
    """Fetches Canonical SMILES from PubChem PUG REST API safely."""
    clean_name = urllib.parse.quote(query_name.strip())
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{clean_name}/property/CanonicalSMILES/JSON"
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    
    try:
        response = requests.get(url, headers=headers, timeout=6)
        if response.status_code == 200:
            props = response.json()["PropertyTable"]["Properties"][0]
            return props.get("CanonicalSMILES")
    except Exception:
        pass
    return None

# --- 5. HEADER SECTION ---
st.markdown("""
<div class="main-header">
    <h1>🧪 FormuAI Computational Engine</h1>
    <p>Real-Time Molecular Property Profiling & Phytochemical Informatics</p>
</div>
""", unsafe_allow_html=True)

# --- 6. SIDEBAR INPUT CONTROLS ---
st.sidebar.header("1. Compound Input Mode")
st.sidebar.markdown("Select Strategy:")

input_mode = st.sidebar.radio(
    "",
    ["Option 1: Global PubChem Search", "Option 2: Interactive SMILES Studio"],
    index=0
)

if input_mode == "Option 1: Global PubChem Search":
    search_query = st.sidebar.text_input("Enter Any Global Drug Name:", value=st.session_state.compound_name)
    if st.sidebar.button("Fetch Compound from PubChem", type="primary", use_container_width=True):
        with st.spinner("Connecting to PubChem..."):
            fetched_smiles = fetch_pubchem_data(search_query)
            if fetched_smiles:
                st.session_state.active_smiles = fetched_smiles
                st.session_state.compound_name = search_query.lower()
                st.sidebar.success(f"Fetched '{search_query}' successfully!")
                st.rerun()
            else:
                st.sidebar.error(f"Could not fetch '{search_query}'. Using current active structure.")

else:
    manual_name = st.sidebar.text_input("Compound Identifier:", value="Custom Structure")
    manual_smiles = st.sidebar.text_area("Paste Canonical SMILES:", value=st.session_state.active_smiles)
    if st.sidebar.button("Apply SMILES", type="primary", use_container_width=True):
        st.session_state.active_smiles = manual_smiles.strip()
        st.session_state.compound_name = manual_name
        st.rerun()

# --- 7. DASHBOARD MAIN DISPLAY ---
st.subheader(f"Active Compound: {st.session_state.compound_name.capitalize()}")

# Active Canonical SMILES Display
current_smiles = st.text_input("Active Canonical SMILES:", value=st.session_state.active_smiles)

if current_smiles != st.session_state.active_smiles:
    st.session_state.active_smiles = current_smiles
    st.rerun()

# --- 8. RDKIT LOCAL CALCULATIONS & VISUALIZATION ---
mol = Chem.MolFromSmiles(st.session_state.active_smiles)

if mol:
    mw = round(Descriptors.MolWt(mol), 2)
    logp = round(Crippen.MolLogP(mol), 2)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
    ro5_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    # Display Property Metrics
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.markdown(f'<div class="metric-card"><div class="metric-val">{mw}</div><div class="metric-lbl">Molecular Weight (g/mol)</div></div>', unsafe_allow_html=True)
    with c2:
        st.markdown(f'<div class="metric-card"><div class="metric-val">{logp}</div><div class="metric-lbl">LogP Partition</div></div>', unsafe_allow_html=True)
    with c3:
        st.markdown(f'<div class="metric-card"><div class="metric-val">{hbd} / {hba}</div><div class="metric-lbl">HBD / HBA</div></div>', unsafe_allow_html=True)
    with c4:
        st.markdown(f'<div class="metric-card"><div class="metric-val">{ro5_violations}</div><div class="metric-lbl">Lipinski Ro5 Violations</div></div>', unsafe_allow_html=True)

    st.write("---")

    # Display 2D Molecular Graph Render
    col_img, col_info = st.columns([1, 1])
    with col_img:
        st.markdown("### 🧬 2D Structure Render")
        img = Draw.MolToImage(mol, size=(450, 450))
        st.image(img, use_column_width=True)

    with col_info:
        st.markdown("### 📋 Pharmacokinetic Profile")
        st.write(f"*TPSA (Polar Surface Area):* {tpsa} Å²")
        
        if ro5_violations == 0:
            st.success("✅ *Lipinski Compliant:* High predicted oral bioavailability.")
        elif ro5_violations == 1:
            st.warning("⚠️ *Moderate Compliance:* 1 Rule violation detected.")
        else:
            st.error(f"❌ *Poor Drug-likeness:* {ro5_violations} Rule violations.")
            
        st.info("💡 Copy the SMILES string above for secondary docking and ADMET evaluations.")
else:
    st.error("Invalid SMILES structure detected. Please double-check structure syntax.")

# --- 9. OWNERSHIP FOOTER ---
st.markdown("""
<div class="footer-container">
    <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap;">
        <div>
            <strong>FormuAI-QbD Engine</strong> | Computational Chemistry & Molecular Profiling
        </div>
        <div>
            <strong>Lead Researcher & Owner:</strong> Mohan Raj Perumal
        </div>
    </div>
</div>
""", unsafe_allow_html=True)
)
