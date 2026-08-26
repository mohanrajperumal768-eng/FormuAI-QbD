import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error

st.set_page_config(page_title="FormuAI-QbD (Advanced)", layout="wide")

# Header Section
st.title("💊 FormuAI-QbD: Advanced Decision-Support Platform")
st.caption("Computer-Aided Pharmaceutical Development & Explainable AI Engine")
st.warning("⚠️ Research Disclaimer: Generated outputs represent virtual candidates for decision support and require empirical laboratory validation.")

# Load Rules & Drug Intelligence Dataset
@st.cache_data
def load_drug_database():
    try:
        return pd.read_csv("drug_rules.csv")
    except Exception:
        return pd.DataFrame({
            "Drug_Name": ["Metformin HCl", "Paracetamol", "Amlodipine Besylate", "Omeprazole"],
            "MW": [165.62, 151.16, 567.1, 345.4],
            "LogP": [-1.43, 0.46, 3.0, 2.23],
            "BCS_Class": ["Class I", "Class III", "Class I", "Class II"],
            "Half_Life_hrs": [6.0, 2.5, 30.0, 1.0],
            "Dose_mg": [500, 500, 5, 20],
            "Solubility": ["Highly Soluble", "Freely Soluble", "Slightly Soluble", "Slightly Soluble"]
        })

df_drugs = load_drug_database()

# Sidebar Navigation
mode = st.sidebar.radio("User Mode", ["Student Mode (Simplified)", "Researcher Mode (Advanced)"])
tabs = st.sidebar.radio("Navigation Menu", [
    "1. Drug Intelligence",
    "2. Dosage Form Ranker",
    "3. Multi-Design Formulation Optimizer",
    "4. Release Kinetics Engine",
    "5. Quantitative QbD Risk Assessment"
])

# TAB 1: DRUG INTELLIGENCE
if tabs == "1. Drug Intelligence":
    st.header("1. API Property & Formulation Significance Analysis")
    selected_drug = st.selectbox("Select API:", df_drugs["Drug_Name"].unique())
    drug_data = df_drugs[df_drugs["Drug_Name"] == selected_drug].iloc[0]

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Physicochemical Profile")
        st.write(f"*Molecular Weight:* {drug_data['MW']} g/mol")
        st.write(f"*LogP:* {drug_data['LogP']}")
        st.write(f"*BCS Classification:* {drug_data['BCS_Class']}")
        st.write(f"*Elimination Half-Life:* {drug_data['Half_Life_hrs']} hrs")
        st.write(f"*Standard Dose:* {drug_data['Dose_mg']} mg")
    
    with col2:
        st.subheader("Formulation Significance")
        if drug_data['Half_Life_hrs'] < 4.0:
            st.info("💡 *Short Half-Life:* High candidate for Extended-Release (ER) matrix design to prevent frequent dosing.")
        if drug_data['BCS_Class'] in ["Class II", "Class IV"]:
            st.warning("⚠️ *Low Solubility:* High risk of dissolution-rate limited bioavailability; consider micronization or lipid carriers.")
        if drug_data['Dose_mg'] >= 500:
            st.error("🚨 *High Dose Load:* Limited suitability for Orally Disintegrating Tablets (ODT) due to tablet mass constraints.")

# TAB 2: DOSAGE FORM RANKER
elif tabs == "2. Dosage Form Ranker":
    st.header("2. Evidence-Based Dosage Form Suitability Ranking")
    
    col1, col2 = st.columns(2)
    with col1:
        target_onset = st.selectbox("Target Action/Release", ["Immediate", "Prolonged (12h)", "Delayed (Gastric Resistance)"])
        swallowing_diff = st.selectbox("Patient Swallowing Difficulty?", ["No", "Yes"])
    with col2:
        dose_val = st.number_input("Unit Dose (mg)", value=250)
        t_half = st.number_input("Drug Half-Life (hrs)", value=3.0)

    if st.button("Evaluate Dosage Forms"):
        rankings = []
        # Evaluation Logic
        er_score = 90 if target_onset == "Prolonged (12h)" and t_half < 6 else 40
        odt_score = 85 if swallowing_diff == "Yes" and dose_val < 200 else 25
        ir_score = 80 if target_onset == "Immediate" else 50
        enteric_score = 88 if target_onset == "Delayed (Gastric Resistance)" else 20

        rankings.append({"Dosage Form": "Extended-Release Matrix Tablet", "Suitability Score": er_score, "Primary Rationale": "Maintains controlled plasma concentrations for short half-life drugs."})
        rankings.append({"Dosage Form": "Orally Disintegrating Tablet (ODT)", "Suitability Score": odt_score, "Primary Rationale": "Enhances compliance for dysphagic patients; requires low unit dose."})
        rankings.append({"Dosage Form": "Immediate-Release Tablet", "Suitability Score": ir_score, "Primary Rationale": "Standard manufacturing, suitable for rapid therapeutic onset."})
        rankings.append({"Dosage Form": "Enteric-Coated Delayed Release", "Suitability Score": enteric_score, "Primary Rationale": "Protects acid-labile APIs from gastric degradation."})

        df_rank = pd.DataFrame(rankings).sort_values(by="Suitability Score", ascending=False)
        st.dataframe(df_rank, use_container_width=True)

# TAB 3: MULTI-DESIGN OPTIMIZER
elif tabs == "3. Multi-Design Formulation Optimizer & XAI":
    st.header("3. Machine Learning Predictor & Virtual Formulation Optimizer")
    
    col1, col2 = st.columns(2)
    with col1:
        target_release = st.slider("Target 12-Hour Cumulative Release (%)", 50, 100, 85)
        polymer_type = st.selectbox("Polymer Selection", ["HPMC K100M", "HPMC K4M", "Eudragit RS PO"])
    with col2:
        compression_force = st.slider("Compression Force (kN)", 5, 25, 12)
        tablet_wt = st.number_input("Target Weight (mg)", value=400)

    # Train Synthetic Model representing literature space
    np.random.seed(42)
    X_train = np.random.uniform(low=[5, 1, 100], high=[40, 10, 500], size=(200, 3))
    # Release dependent on polymer %, binder %, and weight
    y_release = 100 - (X_train[:, 0] * 1.5) - (X_train[:, 1] * 0.5) + np.random.normal(0, 2, 200)
    model = RandomForestRegressor(n_estimators=100, random_state=42).fit(X_train, y_release)

    # Generate Virtual Candidate Designs
    designs = []
    for pol_conc in [10, 18, 25, 32]:
        for binder_conc in [2, 5]:
            pred = model.predict([[pol_conc, binder_conc, tablet_wt]])[0]
            diff = abs(pred - target_release)
            designs.append({
                "Design ID": f"DES-{pol_conc}-{binder_conc}",
                "Polymer %": pol_conc,
                "Binder %": binder_conc,
                "Predicted 12h Release (%)": round(pred, 2),
                "Deviation from Target (%)": round(diff, 2)
            })
    
    df_designs = pd.DataFrame(designs).sort_values(by="Deviation from Target (%)").head(3)
    st.subheader("Top 3 Virtual Design Candidates")
    st.table(df_designs)

    if mode == "Researcher Mode (Advanced)":
        st.subheader("Explainable AI (Feature Importance)")
        importances = model.feature_importances_
        fig_feat = px.bar(x=["Polymer Concentration", "Binder Concentration", "Tablet Weight"], y=importances,
                          labels={'x': 'Formulation Parameter', 'y': 'Relative Feature Importance'}, title="Global Model Feature Drivers")
        st.plotly_chart(fig_feat)

# TAB 4: RELEASE KINETICS ENGINE
elif tabs == "4. Release Kinetics Engine":
    st.header("4. Kinetic Model Fitting & Dissolution Modeling")
    st.write("Input empirical or virtual dissolution profiles to establish release mechanisms.")

    time_str = st.text_input("Time Points (hours, comma separated):", "1, 2, 4, 6, 8, 12")
    rel_str = st.text_input("Cumulative Drug Release (%, comma separated):", "15, 28, 48, 65, 78, 90")

    try:
        t = np.array([float(x.strip()) for x in time_str.split(",")])
        r = np.array([float(x.strip()) for x in rel_str.split(",")])

        # Mathematical Fitting
        r2_zero = r2_score(r, np.poly1d(np.polyfit(t, r, 1))(t))
        r2_first = r2_score(np.log(100 - r + 1e-5), np.poly1d(np.polyfit(t, np.log(100 - r + 1e-5), 1))(t))
        r2_higuchi = r2_score(r, np.poly1d(np.polyfit(np.sqrt(t), r, 1))(np.sqrt(t)))

        # Korsmeyer-Peppas (log t vs log Mt/Minf for release <= 60%)
        mask = r <= 60
        if sum(mask) >= 2:
            n_val, _ = np.polyfit(np.log(t[mask]), np.log(r[mask] / 100 + 1e-5), 1)
        else:
            n_val = 0.45

        # Visual Plotting
        fig_kin = go.Figure()
        fig_kin.add_trace(go.Scatter(x=t, y=r, mode='markers+lines', name='Observed Release Profile'))
        fig_kin.update_layout(title="In-Vitro Dissolution Fit", xaxis_title="Time (hours)", yaxis_title="Cumulative Release (%)")
        st.plotly_chart(fig_kin)

        # Kinetic Results Table
        kin_results = pd.DataFrame({
            "Kinetic Model": ["Zero-Order", "First-Order", "Higuchi (Matrix Diffusion)", "Korsmeyer-Peppas"],
            "Metric / Parameter": [f"R² = {r2_zero:.4f}", f"R² = {r2_first:.4f}", f"R² = {r2_higuchi:.4f}", f"Release Exponent (n) = {n_val:.3f}"],
            "Mechanism Rationale": [
                "Concentration-independent release",
                "Concentration-dependent release",
                "Fickian diffusion through planar matrix",
                "Non-Fickian transport / anomalous diffusion" if n_val > 0.45 else "Fickian diffusion"
            ]
        })
        st.table(kin_results)

    except Exception as e:
        st.error(f"Error executing kinetic analysis: {e}. Please enter matching numeric datasets.")

# TAB 5: QUANTITATIVE QBD RISK ASSESSMENT
elif tabs == "5. Quantitative QbD Risk Assessment":
    st.header("5. Quantitative Quality by Design (QbD) Risk Priority Matrix")
    st.write("Calculate Risk Priority Numbers ($RPN = Severity \times Occurrence \times Detectability$).")

    # Interactive Risk Table
    risk_data = pd.DataFrame({
        "Process/Formulation Parameter": ["Polymer Concentration", "Compression Force", "Mixing / Blending Time", "Particle Size Distribution"],
        "Critical Quality Attribute (CQA)": ["12h Drug Release Profile", "Tablet Friability & Hardness", "Content Uniformity", "Dissolution Rate"],
        "Severity (1-10)": [9, 7, 8, 8],
        "Occurrence (1-10)": [6, 4, 3, 5],
        "Detectability (1-10)": [4, 3, 5, 4]
    })

    edited_df = st.data_editor(risk_data, num_rows="dynamic")
    
    # Calculate RPN
    edited_df["RPN"] = edited_df["Severity (1-10)"] * edited_df["Occurrence (1-10)"] * edited_df["Detectability (1-10)"]
    
    def classify_risk(rpn):
        if rpn >= 200:
            return "High Risk (Requires Control Strategy)"
        elif rpn >= 100:
            return "Medium Risk"
        else:
            return "Low Risk"

    edited_df["Risk Category"] = edited_df["RPN"].apply(classify_risk)
    
    st.subheader("Evaluated Risk Matrix")
    st.dataframe(edited_df.sort_values(by="RPN", ascending=False), use_container_width=True)

    # Export Summary
    csv = edited_df.to_csv(index=False)
    st.download_button("Download Executive QbD Audit Report (CSV)", csv, "FormuAI_QbD_Audit_Report.csv", "text/csv")
