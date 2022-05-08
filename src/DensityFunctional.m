(* ::Section::Closed:: *)
(*KineticEnergyFunctional*)


(* expression for kinetic energy functional *)
KineticEnergyFunctional[density_, potential_]:=Module[
    {kineticEnergy},

    (* define the expression *)
    kineticEnergy = 1./2 * Dot[r, Grad] * potential * density;

    (* return the integral *)
    Integrate[kineticEnergy, r];
];


(*::Section::Closed::*)
(*PotentialEnergyFunctional*)


(* define the simple form of the potential energy functional *)
PotentialEnergyFunctional[density_, structure_]:=Module[
    {hartreePotential, xcPotential, nuclearPotential},

    (* find the hartree potential *)
    hartreePotential = HartreeFunctional[density];

    (* find the exchange-correlation potential *)
    xcPotential = XCFunctional[density];

    (* find the potential due to the Coulomb attraction with nuclei *)
    nuclearPotential = NuclearFunctional[density, structure];

    (* return the integral of the sum *)
    Integrate[(hartreePotential + xcPotential + nuclearPotential) * density[r], r]
];


(* hartree functional *)
HartreeFunctional[density_]:=Module[
    {},

    (* define the integrand *)
    Integrate[density[rPrime]/Norm[r-rPrime], rPrime];
];


(* exchange-correlation functional *)
XCFunctional[density_]:=XCFunctional[density, "PBE"];
XCFunctional[density_, type_String]:=Module[
    {densityFile, outputCode},

    (* convert density to file *)
    densityFile = createDensityFile[density];

    (* switch behavior based on type *)
    outputCode = Switch[type,
        "PBE", Run["python pbe_libxc.py -d "<>densityFile],
        _, Message[Error::UnknownXCFunctional, type]; $Failed
    ];

    (* read output file and convert to a function *)
    potentialData = readPotentialData[$LibXCOutputFileName];
    CreateFunction[potentialData]
];


(* nuclear potential functional *)
NuclearFunctional[density_, structure_]:=Module[
    {charges, positions, potentialContributions},

    (* unpack the charges and atomic positions from structure object *)
    (* structure has the form of {<|Z->charge, Position->{x,y,z}|>,...}*)
    charges = Lookup[structure, Z];
    positions = Lookup[structure, Position];

    (* Map thread these values into correct Coulomb potential energy form *)
    potentialContributions = MapThread[-#1/Norm[r-#2]&, {charges, positions}];

    (* return the sum *)
    Total[potentialContributions]
];


(* Create function from data using interpolate *)
CreateFunction[data_]:=Module[
    {},

    (* use interpolate to create an interpolating function of the data *)
    (* data is in the format of {{r, value}, ... }*)
    (* as a result, it can be used directly with the Interpolation function *)
    Interpolation[data]
];


ReadData[fileName_String]:=Module[
    {},

    (* import the file *)
    (* impose HDF5 format to ensure that large datasets are handled properly *)
    Import[fileName, "HDF5"];
];
