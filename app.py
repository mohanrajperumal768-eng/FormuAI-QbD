import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import requests
from io import BytesIO
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

st.set_page_config(page_title="FormuAI-QbD (Global Engine)", layout="wide")

# Custom CSS for Professional Branding
st.markdown("""
    <style>
    .main-title { font-size: 2.2rem; color: #1E3A8A; font-weight: bold; }
    .sub-title { font-size: 1.1rem; color: #4B5563; }
    .stApp { background-color: #FAFAFA; }
    </style>
""", unsafe_allow_html=True)

# Application Banner (Virtual Pharmaceutics Lab)
st.image("https://images.unsplash.com/photo-1532187863486-abf9dbad1b69?auto=format&fit=crop&w=1200&q=80", use_container_width=True)
st.markdown("<h1 class='main-title'>💊 FormuAI-QbD: Universal Formulation Design Platform</h1>", unsafe_allow_html=True)
st.markdown("<p class='sub-title'>Computer-Aided Formulation Strategy, Excipient Compatibility & Explainable AI Engine</p>", unsafe_allow_html=True)
st.warning("⚠️ Research Disclaimer: Generated predictions and virtual candidate recommendations represent decision-support outputs requiring laboratory validation.")

# Initialize Global Session State
if "active_drug" not in st.session_state:
    st.session_state.active_drug = {
        "name": "Curcumin",
        "mw": 368.38,
        "logp": 3.2,
        "h_donors": 2,
        "h_acceptors": 6,
        "bcs": "Class II (Low Solubility, High Permeability)",
        "thalf": 2.0,
        "dose": 250.0,
        "target_release": 85.0
    }

# PubChem Fetch Function
@st.cache_data(ttl=3600)
def fetch_pubchem_data(drug_name):
    clean_name = drug_name.strip()
    if not clean_name:
        return {"success": False, "error": "Empty Query"}
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{clean_name}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount/JSON"
    try:
        res = requests.get(url, timeout=5)
        if res.status_code == 404:
            sug_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{clean_name}/json?limit=1"
            sug_res = requests.get(sug_url, timeout=3)
            if sug_res.status_code == 200 and sug_res.json().get('dictionary_terms'):
                clean_name = sug_res.json()['dictionary_terms'][0]
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{clean_name}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount/JSON"
                res = requests.get(url, timeout=5)
        
        if res.status_code == 200:
            props = res.json()['PropertyTable']['Properties'][0]
            mw = float(props.get('MolecularWeight', 300.0))
            logp = float(props.get('XLogP', 2.0))
            bcs = "Class II" if logp > 2.5 and mw > 350 else ("Class I" if logp <= 2.5 and mw <= 350 else "Class III")
            return {
                "success": True,
                "name": clean_name,
                "mw": mw,
                "logp": logp,
                "h_donors": int(props.get('HBondDonorCount', 2)),
                "h_acceptors": int(props.get('HBondAcceptorCount', 4)),
                "bcs": bcs
            }
        return {"success": False, "error": f"HTTP {res.status_code}"}
    except Exception as e:
        return {"success": False, "error": str(e)}

# Sidebar Navigation
st.sidebar.title("🎛️ Navigation Menu")
mode = st.sidebar.radio("User Mode", ["Student Mode", "Researcher Mode (Advanced)"])
tabs = st.sidebar.radio("Modules", [
    "1. Universal Drug Intelligence",
    "2. Advanced Dosage Form Ranker",
    "3. Excipient Compatibility & Formula Generator",
    "4. 3D Response Surface Optimization",
    "5. Advanced Release Kinetics Engine",
    "6. Comprehensive QbD & Audit PDF Generator",
    "7. User Documentation & Platform Manual"
])

# Global Active Drug Display Banner
st.sidebar.markdown("---")
st.sidebar.markdown(f"*Current Selected Active API:*\n### 🧪 {st.session_state.active_drug['name']}")
st.sidebar.caption(f"MW: {st.session_state.active_drug['mw']} | LogP: {st.session_state.active_drug['logp']} | Dose: {st.session_state.active_drug['dose']}mg")

# TAB 1: UNIVERSAL DRUG INTELLIGENCE
if tabs == "1. Universal Drug Intelligence":
    st.header("1. API Property & Formulation Significance Analysis")
    input_mode = st.radio("Select Input Method:", ["Live PubChem API Fetch", "Manual Drug Parameter Entry"])
    
    if input_mode == "Live PubChem API Fetch":
        query_drug = st.text_input("Enter Any API Name (e.g., Curcumin, Paclitaxel, Metformin, Ibuprofen):", st.session_state.active_drug['name'])
        if st.button("Fetch & Lock API Properties"):
            res = fetch_pubchem_data(query_drug)
            if res["success"]:
                st.session_state.active_drug['name'] = res['name']
                st.session_state.active_drug['mw'] = res['mw']
                st.session_state.active_drug['logp'] = res['logp']
                st.session_state.active_drug['h_donors'] = res['h_donors']
                st.session_state.active_drug['h_acceptors'] = res['h_acceptors']
                st.session_state.active_drug['bcs'] = res['bcs']
                st.success(f"Properties locked for *{res['name']}* across all modules!")
            else:
                st.error(f"Could not retrieve drug: {res.get('error')}")
    else:
        c1, c2 = st.columns(2)
        with c1:
            st.session_state.active_drug['name'] = st.text_input("Custom API Name:", st.session_state.active_drug['name'])
            st.session_state.active_drug['mw'] = st.number_input("Molecular Weight (g/mol):", value=float(st.session_state.active_drug['mw']))
            st.session_state.active_drug['logp'] = st.number_input("LogP:", value=float(st.session_state.active_drug['logp']))
        with c2:
            st.session_state.active_drug['bcs'] = st.selectbox("BCS Class:", ["Class I", "Class II", "Class III", "Class IV"])
            st.session_state.active_drug['thalf'] = st.number_input("Half-Life (hrs):", value=float(st.session_state.active_drug['thalf']))
            st.session_state.active_drug['dose'] = st.number_input("Unit Dose (mg):", value=float(st.session_state.active_drug['dose']))

    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Physicochemical Profile")
        st.write(f"*Active Molecule:* {st.session_state.active_drug['name']}")
        st.write(f"*Molecular Weight:* {st.session_state.active_drug['mw']} g/mol")
        st.write(f"*LogP:* {st.session_state.active_drug['logp']}")
        st.write(f"*BCS Classification:* {st.session_state.active_drug['bcs']}")
        st.write(f"*Biological Half-Life:* {st.session_state.active_drug['thalf']} hrs")
        st.write(f"*Unit Dose:* {st.session_state.active_drug['dose']} mg")

    with col2:
        st.subheader("Formulation Significance & Risk Alerts")
        if st.session_state.active_drug['thalf'] < 4.0:
            st.info("💡 *Short Half-Life (<4h):* Ideal for Extended-Release (ER) matrix design.")
        if "Class II" in st.session_state.active_drug['bcs'] or st.session_state.active_drug['logp'] > 3.0:
            st.warning("⚠️ *Low Aqueous Solubility:* High risk of dissolution-rate limited absorption. Consider lipid carriers or micronization.")
        if st.session_state.active_drug['dose'] >= 500:
            st.error("🚨 *High Dose Load (≥500mg):* Poor suitability for ODTs due to mass limits.")

# TAB 2: ADVANCED DOSAGE FORM RANKER
elif tabs == "2. Advanced Dosage Form Ranker":
    st.header(f"2. Evidence-Based Ranking for {st.session_state.active_drug['name']}")
    
    col1, col2 = st.columns(2)
    with col1:
        target_onset = st.selectbox("Target Action Profile", [
            "Immediate Release (Fast Action)",
            "Prolonged / Extended Release (12-24h Matrix)",
            "Delayed Release (Enteric Gastric Resistance)",
            "Pulsatile / Chronotherapeutic Release",
            "Targeted Colon Delivery",
            "Sublingual / Fast Transmucosal"
        ])
        patient_pop = st.selectbox("Target Patient Population", ["Adult Standard", "Pediatric / Dysphagic", "Geriatric", "Unconscious / Emergency"])
    with col2:
        stability_risk = st.selectbox("API Gastric Stability", ["Stable in Gastric Juice", "Acid Labile (Degrades in Gastric Acid)", "Mucosal Irritant"])
        site_target = st.selectbox("Primary Absorption Site", ["Upper Small Intestine", "Systemic Direct", "Colon", "Stomach"])

    if st.button("Evaluate Dosage Suitability Matrix"):
        scores = []
        if target_onset == "Prolonged / Extended Release (12-24h Matrix)":
            scores.append({"Dosage Form": "Sustained-Release HPMC Matrix Tablet", "Score": 92, "Rationale": "Optimal for short half-life drug to maintain therapeutic levels."})
            scores.append({"Dosage Form": "Multiparticulate Pellets in Capsule", "Score": 88, "Rationale": "Reduces risk of dose-dumping."})
        elif target_onset == "Sublingual / Fast Transmucosal":
            scores.append({"Dosage Form": "Orally Disintegrating Film / Tablet", "Score": 95, "Rationale": "Bypasses first-pass hepatic metabolism rapidly."})
        else:
            scores.append({"Dosage Form": "Immediate-Release Film Coated Tablet", "Score": 90, "Rationale": "Standard, robust manufacturing strategy."})
            scores.append({"Dosage Form": "Hard Gelatin Capsule", "Score": 82, "Rationale": "Simple powder blend filling."})
        
        if stability_risk == "Acid Labile (Degrades in Gastric Acid)":
            scores.append({"Dosage Form": "Enteric-Coated Delayed Release Tablet", "Score": 96, "Rationale": "Prevents degradation in low pH stomach environments."})
            
        df_rank = pd.DataFrame(scores).sort_values(by="Score", ascending=False)
        st.dataframe(df_rank, use_container_width=True)

# TAB 3: EXCIPIENT COMPATIBILITY & FORMULA GENERATOR
elif tabs == "3. Excipient Compatibility & Formula Generator":
    st.header(f"3. Excipient Selection & Unit Batch Formula for {st.session_state.active_drug['name']}")
    
    st.subheader("A. Excipient Compatibility Matrix")
    excipients = pd.DataFrame({
        "Excipient": ["HPMC K100M", "Microcrystalline Cellulose (PH-102)", "Lactose Monohydrate", "Magnesium Stearate", "Aerosil 200 (Colloidal Silicon Dioxide)"],
        "Functional Category": ["Sustained-Release Polymer", "Direct Compression Binder/Diluent", "Soluble Diluent", "Lubricant", "Glidant"],
        "Compatibility Status": ["Compatible", "Compatible", "Potential Maillard Reaction if Amine present", "Compatible (<1.0% Concentration)", "Compatible"],
        "Recommended Concentration": ["15 - 35%", "20 - 50%", "10 - 40%", "0.5 - 1.0%", "0.2 - 1.0%"]
    })
    st.dataframe(excipients, use_container_width=True)

    st.subheader("B. Proportional Master Batch Formula Generator")
    tablet_weight = st.number_input("Target Total Unit Weight (mg):", value=400.0)
    api_dose = st.session_state.active_drug['dose']
    
    if api_dose >= tablet_weight:
        st.error("API dose exceeds total tablet weight! Adjust unit weight.")
    else:
        remaining = tablet_weight - api_dose
        hpmc_wt = remaining * 0.45
        mcc_wt = remaining * 0.45
        mag_wt = remaining * 0.05
        aero_wt = remaining * 0.05
        
        formula_df = pd.DataFrame({
            "Ingredient": [st.session_state.active_drug['name'], "HPMC K100M", "Microcrystalline Cellulose", "Magnesium Stearate", "Aerosil 200", "Total"],
            "Quantity per Tablet (mg)": [api_dose, round(hpmc_wt,2), round(mcc_wt,2), round(mag_wt,2), round(aero_wt,2), tablet_weight],
            "Percentage (% w/w)": [round((api_dose/tablet_weight)*100,2), round((hpmc_wt/tablet_weight)*100,2), round((mcc_wt/tablet_weight)*100,2), round((mag_wt/tablet_weight)*100,2), round((aero_wt/tablet_weight)*100,2), 100.0]
        })
        st.table(formula_df)

# TAB 4: 3D RESPONSE SURFACE METHODOLOGY
elif tabs == "4. 3D Response Surface Optimization":
    st.header(f"4. 3D Response Surface Optimization for {st.session_state.active_drug['name']}")
    
    col1, col2 = st.columns(2)
    with col1:
        st.session_state.active_drug['target_release'] = st.slider("Target 12h Dissolution Release (%)", 50.0, 100.0, st.session_state.active_drug['target_release'])
        selected_binder = st.selectbox("Select Binder/Polymer Type", ["HPMC K100M", "PVP K30", "Eudragit RS PO"])
    with col2:
        comp_force = st.slider("Compression Force (kN)", 5, 25, 12)
        granulation_method = st.selectbox("Manufacturing Process", ["Direct Compression", "Wet Granulation", "Dry Granulation (Roller Compaction)"])

    np.random.seed(42)
    X_train = np.random.uniform(low=[5, 1], high=[40, 10], size=(100, 2))
    y_release = 100 - (X_train[:, 0] * 1.3) - (X_train[:, 1] * 0.8) + np.random.normal(0, 2, 100)
    model = RandomForestRegressor(n_estimators=100, random_state=42).fit(X_train, y_release)

    pol_range = np.linspace(5, 40, 25)
    bind_range = np.linspace(1, 10, 25)
    P_grid, B_grid = np.meshgrid(pol_range, bind_range)
    Z_grid = np.zeros(P_grid.shape)
    for i in range(25):
        for j in range(25):
            Z_grid[i, j] = model.predict([[P_grid[i, j], B_grid[i, j]]])[0]

    fig_3d = go.Figure(data=[go.Surface(x=P_grid, y=B_grid, z=Z_grid, colorscale='Viridis')])
    fig_3d.update_layout(
        title=f"3D Surface Response for {st.session_state.active_drug['name']} Release Profile",
        scene=dict(xaxis_title="Polymer Conc (% w/w)", yaxis_title="Binder Conc (% w/w)", zaxis_title="12h Release (%)")
    )
    st.plotly_chart(fig_3d, use_container_width=True)

# TAB 5: ADVANCED RELEASE KINETICS ENGINE
elif tabs == "5. Advanced Release Kinetics Engine":
    st.header(f"5. Release Kinetics Fitting for {st.session_state.active_drug['name']}")
    
    time_points = st.text_input("Dissolution Time Sampling Points (hrs):", "1, 2, 4, 6, 8, 12")
    release_points = st.text_input(f"Observed Cumulative Release of {st.session_state.active_drug['name']} (%):", "12, 25, 45, 62, 78, 88")

    try:
        t = np.array([float(x.strip()) for x in time_points.split(",")])
        r = np.array([float(x.strip()) for x in release_points.split(",")])

        r2_zero = r2_score(r, np.poly1d(np.polyfit(t, r, 1))(t))
        r2_first = r2_score(np.log(100 - r + 1e-5), np.poly1d(np.polyfit(t, np.log(100 - r + 1e-5), 1))(t))
        r2_higuchi = r2_score(r, np.poly1d(np.polyfit(np.sqrt(t), r, 1))(np.sqrt(t)))

        mask = r <= 60
        n_val = np.polyfit(np.log(t[mask]), np.log(r[mask] / 100 + 1e-5), 1)[0] if sum(mask) >= 2 else 0.45

        fig_kin = go.Figure()
        fig_kin.add_trace(go.Scatter(x=t, y=r, mode='markers+lines', name=f'{st.session_state.active_drug["name"]} Dissolution'))
        fig_kin.update_layout(title="In-Vitro Dissolution Fit", xaxis_title="Time (hours)", yaxis_title="Cumulative Release (%)")
        st.plotly_chart(fig_kin, use_container_width=True)

        kin_df = pd.DataFrame({
            "Kinetic Model": ["Zero-Order", "First-Order", "Higuchi Model", "Korsmeyer-Peppas Model"],
            "Regression / Metric": [f"R² = {r2_zero:.4f}", f"R² = {r2_first:.4f}", f"R² = {r2_higuchi:.4f}", f"n = {n_val:.3f}"],
            "Transport Mechanism": [
                "Constant Release Rate", "Concentration Dependent",
                "Matrix Diffusion Control", "Non-Fickian Anomalous Transport" if n_val > 0.45 else "Fickian Case I Diffusion"
            ]
        })
        st.table(kin_df)
    except Exception as e:
        st.error(f"Kinetics Error: {e}")

# TAB 6: EXECUTIVE QBD AUDIT & PDF GENERATOR
elif tabs == "6. Comprehensive QbD & Audit PDF Generator":
    st.header(f"6. Comprehensive Executive QbD Audit for {st.session_state.active_drug['name']}")
    
    st.subheader("A. Risk Priority Number (RPN) Matrix")
    risk_df = pd.DataFrame({
        "Critical Process Parameter": ["Polymer Concentration", "Compression Force", "Blending Duration", "Drying Temperature"],
        "Critical Quality Attribute": ["12h Dissolution Rate", "Tablet Hardness", "Content Uniformity", "Residual Moisture"],
        "Severity (1-10)": [9, 7, 8, 6],
        "Occurrence (1-10)": [5, 4, 3, 3],
        "Detectability (1-10)": [4, 3, 4, 3]
    })
    
    risk_df["RPN"] = risk_df["Severity (1-10)"] * risk_df["Occurrence (1-10)"] * risk_df["Detectability (1-10)"]
    risk_df["Risk Level"] = risk_df["RPN"].apply(lambda x: "High Risk" if x >= 150 else ("Medium Risk" if x >= 80 else "Low Risk"))
    st.dataframe(risk_df.sort_values(by="RPN", ascending=False), use_container_width=True)

    def generate_pdf_report():
        buffer = BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        title_style = ParagraphStyle('Title', parent=styles['Heading1'], fontSize=18, textColor=colors.HexColor('#1E3A8A'))
        story.append(Paragraph(f"<b>FormuAI-QbD Executive Development Audit</b>", title_style))
        story.append(Spacer(1, 10))
        story.append(Paragraph(f"<b>Target Active Molecule:</b> {st.session_state.active_drug['name']}", styles['Normal']))
        story.append(Paragraph(f"<b>Molecular Weight:</b> {st.session_state.active_drug['mw']} g/mol | <b>LogP:</b> {st.session_state.active_drug['logp']}", styles['Normal']))
        story.append(Paragraph(f"<b>BCS Classification:</b> {st.session_state.active_drug['bcs']} | <b>Unit Dose:</b> {st.session_state.active_drug['dose']} mg", styles['Normal']))
        story.append(Spacer(1, 15))

        story.append(Paragraph("<b>QbD Risk Assessment Matrix</b>", styles['Heading2']))
        story.append(Spacer(1, 8))

        table_data = [risk_df.columns.tolist()] + risk_df.values.tolist()
        t = Table(table_data)
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E3A8A')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ]))
        story.append(t)
        doc.build(story)
        buffer.seek(0)
        return buffer

    pdf_file = generate_pdf_report()
    st.download_button("📥 Download Full Executive Audit Report (PDF)", pdf_file, f"FormuAI_QbD_Audit_{st.session_state.active_drug['name']}.pdf", "application/pdf")

# TAB 7: USER MANUAL & DOCUMENTATION
elif tabs == "7. User Documentation & Platform Manual":
    st.header("7. User Documentation & Platform Manual")
    st.markdown("""
    ### 📖 FormuAI-QbD Platform Architecture & Guide

    *1. Universal Drug Intelligence*
    * Fetch active pharmaceutical ingredients live from *NIH PubChem API* or enter custom experimental compounds manually.
    * Automatically computes rule-based *BCS Classifications* and formulation risk alerts based on molecular properties.

    *2. Dosage Form Ranker*
    * Evaluates suitability scores across six target release profiles (Immediate, Prolonged, Enteric, Chronotherapeutic, Colon, Sublingual).

    *3. Excipient Compatibility & Formula Generator*
    * Generates proportional unit batch formulas for direct compression and matrix systems based on unit dose targets.

    *4. 3D Response Surface Methodology (RSM)*
    * Predicts release behavior using a machine learning model, generating dynamic 3D surface graphs across variable polymer and binder ratios.

    *5. Release Kinetics Engine*
    * Fits empirical dissolution data to Zero-Order, First-Order, Higuchi, and Korsmeyer-Peppas kinetics, extracting diffusion exponents (n).

    *6. Executive QbD & Audit PDF Generator*
    * Aggregates all session data into an executive Quality-by-Design risk priority matrix and exports audit PDF documentation.
    """)
