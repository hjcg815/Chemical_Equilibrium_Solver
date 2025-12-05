# reactions.py
# List of 20 chemical reactions for equilibrium solver

reactions = [
    {
        "id": 1,
        "name": "Haber process (ammonia synthesis)",
        "stoichiometry": {"N2": -1, "H2": -3, "NH3": 2},
        "equation": "N2 + 3 H2 ⇌ 2 NH3",
        "description": "Used industrially to produce ammonia for fertilizers via nitrogen and hydrogen reaction."
    },
    {
        "id": 2,
        "name": "Water–gas shift",
        "stoichiometry": {"CO": -1, "H2O": -1, "CO2": 1, "H2": 1},
        "equation": "CO + H2O ⇌ CO2 + H2",
        "description": "Increases hydrogen yield from syngas for fuel and chemical production."
    },
    {
        "id": 3,
        "name": "Steam reforming of methane",
        "stoichiometry": {"CH4": -1, "H2O": -1, "CO": 1, "H2": 3},
        "equation": "CH4 + H2O ⇌ CO + 3 H2",
        "description": "Primary industrial method for producing hydrogen and synthesis gas from methane."
    },
    {
        "id": 4,
        "name": "Boudouard reaction",
        "stoichiometry": {"CO": -2, "CO2": 1, "C": 1},
        "equation": "2 CO ⇌ CO2 + C",
        "description": "Important in metallurgical processes; leads to carbon deposition in furnaces."
    },
    {
        "id": 5,
        "name": "Hydrogen–iodine formation",
        "stoichiometry": {"H2": -1, "I2": -1, "HI": 2},
        "equation": "H2 + I2 ⇌ 2 HI",
        "description": "Example reaction for studying equilibrium kinetics and thermodynamics."
    },
    {
        "id": 6,
        "name": "Dinitrogen tetroxide dissociation",
        "stoichiometry": {"N2O4": -1, "NO2": 2},
        "equation": "N2O4 ⇌ 2 NO2",
        "description": "Used in chemical engineering studies to demonstrate gas-phase equilibrium."
    },
    {
        "id": 7,
        "name": "Sulfur dioxide oxidation",
        "stoichiometry": {"SO2": -2, "O2": -1, "SO3": 2},
        "equation": "2 SO2 + O2 ⇌ 2 SO3",
        "description": "Industrial process for sulfuric acid production (Contact Process)."
    },
    {
        "id": 8,
        "name": "Hydrogen–bromine formation",
        "stoichiometry": {"H2": -1, "Br2": -1, "HBr": 2},
        "equation": "H2 + Br2 ⇌ 2 HBr",
        "description": "Demonstrates halogen-hydrogen equilibria; used in lab experiments."
    },
    {
        "id": 9,
        "name": "Hydrogen combustion",
        "stoichiometry": {"H2": -2, "O2": -1, "H2O": 2},
        "equation": "2 H2 + O2 ⇌ 2 H2O",
        "description": "Basic combustion reaction; fundamental in thermodynamics studies."
    },
    {
        "id": 10,
        "name": "Methanation reaction",
        "stoichiometry": {"CO2": -1, "H2": -4, "CH4": 1, "H2O": 2},
        "equation": "CO2 + 4 H2 ⇌ CH4 + 2 H2O",
        "description": "Sabatier process; produces methane from CO2 and H2, important in energy storage."
    },
    {
        "id": 11,
        "name": "Nitric oxide formation",
        "stoichiometry": {"N2": -1, "O2": -1, "NO": 2},
        "equation": "N2 + O2 ⇌ 2 NO",
        "description": "Relevant to high-temperature combustion and NOx formation in engines."
    },
    {
        "id": 12,
        "name": "Methane cracking",
        "stoichiometry": {"CH4": -1, "C": 1, "H2": 2},
        "equation": "CH4 ⇌ C + 2 H2",
        "description": "Used industrially for hydrogen production and carbon deposition processes."
    },
    {
        "id": 13,
        "name": "Ammonia decomposition",
        "stoichiometry": {"NH3": -2, "N2": 1, "H2": 3},
        "equation": "2 NH3 ⇌ N2 + 3 H2",
        "description": "Relevant in ammonia recycling and hydrogen production studies."
    },
    {
        "id": 14,
        "name": "Methanol synthesis",
        "stoichiometry": {"CO": -1, "H2": -2, "CH3OH": 1},
        "equation": "CO + 2 H2 ⇌ CH3OH",
        "description": "Industrial methanol production from syngas for fuels and chemicals."
    },
    {
        "id": 15,
        "name": "Reverse water–gas shift",
        "stoichiometry": {"CO2": -1, "H2": -1, "CO": 1, "H2O": 1},
        "equation": "CO2 + H2 ⇌ CO + H2O",
        "description": "Used to adjust CO/H2 ratio in syngas; produces CO from CO2."
    },
    {
        "id": 16,
        "name": "Hydrogen chloride decomposition",
        "stoichiometry": {"HCl": -2, "H2": 1, "Cl2": 1},
        "equation": "2 HCl ⇌ H2 + Cl2",
        "description": "Industrial relevance in chlor-alkali processes; demonstrates gas-phase equilibrium."
    },
    {
        "id": 17,
        "name": "Carbon disulfide combustion",
        "stoichiometry": {"CS2": -1, "O2": -3, "CO2": 1, "SO2": 2},
        "equation": "CS2 + 3 O2 ⇌ CO2 + 2 SO2",
        "description": "Demonstrates combustion of sulfur-containing compounds in industrial processes."
    },
    {
        "id": 18,
        "name": "Ammonium chloride formation",
        "stoichiometry": {"NH3": -1, "HCl": -1, "NH4Cl": 1},
        "equation": "NH3 + HCl ⇌ NH4Cl",
        "description": "Used in laboratory synthesis and chemical education for solid-gas equilibria."
    },
    {
        "id": 19,
        "name": "Complete CO oxidation",
        "stoichiometry": {"CO": -2, "O2": -1, "CO2": 2},
        "equation": "2 CO + O2 ⇌ 2 CO2",
        "description": "Fundamental reaction in combustion and environmental control studies."
    },
    {
        "id": 20,
        "name": "Ethylene hydrogenation",
        "stoichiometry": {"C2H4": -1, "H2": -1, "C2H6": 1},
        "equation": "C2H4 + H2 ⇌ C2H6",
        "description": "Important in petrochemical industry; converts unsaturated hydrocarbons to saturated ones."
    }
]
