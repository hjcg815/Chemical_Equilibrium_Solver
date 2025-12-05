# app.py
import streamlit as st
from solver import solve_equilibrium
from reactions import reactions

st.set_page_config(page_title="Chemical Equilibrium Solver", layout="wide")
st.title("Chemical Equilibrium Solver")
st.markdown("""
This app calculates equilibrium compositions, ΔH°, ΔS°, ΔG°, equilibrium constant K,
and extent of reaction ξ for selected chemical reactions.
""")

# Select reaction
reaction_names = [r["name"] for r in reactions]
selected_name = st.selectbox("Select a reaction:", reaction_names)

reaction = next((r for r in reactions if r["name"] == selected_name), None)

if reaction:
    # Display balanced chemical equation
    st.subheader("Balanced Chemical Equation")
    st.markdown(f"**{reaction['equation']}**")

    # Display reaction description
    if "description" in reaction:
        st.subheader("Reaction Description / Industrial Relevance")
        st.write(reaction["description"])

    # Input temperature and pressure
    st.subheader("Enter reaction conditions:")
    T = st.number_input("Temperature (K):", min_value=0.0, value=298.15)
    P = st.number_input("Pressure (atm):", min_value=0.0, value=1.0)

    # Input initial moles
    st.subheader("Initial moles of each species:")
    n0 = {}
    for species in reaction["stoichiometry"]:
        n0[species] = st.number_input(f"{species} (mol):", min_value=0.0, value=1.0)

    # Solve equilibrium
    if st.button("Calculate Equilibrium"):
        results = solve_equilibrium(reaction, n0, T, P)

        st.subheader("Thermodynamic Properties:")
        st.write(f"ΔH° = {results['ΔH']:.2f} kJ/mol")
        st.write(f"ΔS° = {results['ΔS']:.2f} kJ/mol·K")
        st.write(f"ΔG° = {results['ΔG']:.2f} kJ/mol")
        st.write(f"Equilibrium constant K = {results['K']:.4e}")
        st.write(f"Extent of reaction ξ = {results['ξ_eq']:.4f} mol")

        st.subheader("Equilibrium Composition:")
        st.table({
            "Species": list(results["n_eq"].keys()),
            "Moles (n_eq)": [f"{v:.4f}" for v in results["n_eq"].values()],
            "Mole Fraction (y_eq)": [f"{v:.4f}" for v in results["y_eq"].values()]
        })
else:
    st.error("Selected reaction not found in reactions.py!")

