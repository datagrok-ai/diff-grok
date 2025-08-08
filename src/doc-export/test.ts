/* eslint-disable max-len */
import {getIVP, IVP} from '../scripting-tools';
import {DER_FORM, FORMULA_TAG} from './format-defs';
import {ModelToDocExporter, ModelToLaTeXExporter} from './model-to-doc-exporter';

const poll = `#name: Pollution
#description: The chemical reaction part of the air pollution model developed at The Dutch National Institute of Public Health and Environmental Protection
#comment: 
  Source: https://archimede.uniba.it/~testset/report/pollu.pdf
#equations:
  dy1/dt = -(r1 + r10 + r14 + r23 + r24) + (r2 + r3 + r9 + r11 + r12 + r22 + r25)

  dy2/dt = -r2 - r3 - r9 - r12 + r1 + r21

  dy3/dt = -r15 + r1 + r17 + r19 + r22

  dy4/dt = -r2 - r16 - r17 - r23 + r15

  dy5/dt = -r3 +2 * r4 + r6 +r7 +r13 + r20

  dy6/dt = -r6 - r8 - r14 - r20 + r3 + 2 * r18

  dy7/dt = -r4 - r5 - r6 + r13

  dy8/dt = r4 + r5 + r6 + r7

  dy9/dt = -r7 - r8

  dy10/dt = -r12 +r7 + r9

  dy11/dt = -r9 - r10 + r8 + r11

  dy12/dt = r9

  dy13/dt = -r11 + r10

  dy14/dt = -r13 + r12

  dy15/dt = r14

  dy16/dt = -r18 - r19 + r16

  dy17/dt = -r20

  dy18/dt = r20

  dy19/dt = -r21 - r22 - r24 + r23 + r25

  dy20/dt = -r25 + r24

#expressions:
  r1 = k1 * y1

  r2 = k2 * y2 * y4

  r3 = k3 * y5 * y2

  r4 = k4 * y7

  r5 = k5 * y7

  r6 = k6 * y7 * y6

  r7 = k7 * y9

  r8 = k8 * y9 * y6

  r9 = k9 * y11 * y2

  r10 = k10 * y11 * y1

  r11 = k11 * y13

  r12 = k12 * y10 * y2

  r13 = k13 * y14

  r14 = k14 * y1 * y6

  r15 = k15 * y3

  r16 = k16 * y4

  r17 = k17 * y4

  r18 = k18 * y16

  r19 = k19 * y16

  r20 = k20 * y17 * y6

  r21 = k21 * y19

  r22 = k22 * y19

  r23 = k23 * y1 * y4

  r24 = k24 * y19 * y1

  r25 = k25 * y20

#argument: t
  t0 = 0   {units: min; caption: Initial; category: Time; min: 0; max: 0.9} [Initial time of simulation]
  t1 = 60  {units: min; caption: Final; category: Time; min: 1; max: 100; step: 1} [Final time of simulation]
  h = 0.1  {units: min; caption: Step; category: Time; min: 0.001; max: 0.1; step: 0.001} [Time step of simulation]


#inits:  
  y1 = 0    {caption: NO2; category: Initial concentrations; min: 0; max: 0.1} [Initial concentration of NO2]
  y2 = 0.2  {caption: NO; category: Initial concentrations; min: 0; max: 0.4} [Initial concentration of NO]
  y3 = 0    {caption: O3P; category: Initial concentrations; min: 0; max: 0.1} [Initial concentration of O3P]
  y4 = 0.04 {caption: O3; category: Initial concentrations; min: 0; max: 0.1} [Initial concentration of O3]
  y5 = 0    {caption: HO2; category: Initial concentrations} [Initial concentration of HO2]
  y6 = 0    {caption: OH; category: Initial concentrations} [Initial concentration of OH]
  y7 = 0.1  {caption: HCHO; category: Initial concentrations} [Initial concentration of HCHO]
  y8 = 0.3  {caption: CO; category: Initial concentrations} [Initial concentration of CO]
  y9 = 0.01 {caption: ALD; category: Initial concentrations} [Initial concentration of ALD] 
  y10 = 0   {caption: MEO2; category: Initial concentrations} [Initial concentration of MEO2]
  y11 = 0   {caption: C2O3; category: Initial concentrations} [Initial concentration of C2O3]
  y12 = 0   {caption: CO2; category: Initial concentrations} [Initial concentration of CO2]
  y13 = 0   {caption: PAN; category: Initial concentrations} [Initial concentration of PAN]
  y14 = 0   {caption: CH3O; category: Initial concentrations} [Initial concentration of CH3O]
  y15 = 0   {caption: HNO3; category: Initial concentrations} [Initial concentration of HNO3]
  y16 = 0   {caption: O1D; category: Initial concentrations} [Initial concentration of O1D]
  y17 = 0.007 {caption: SO2; category: Initial concentrations} [Initial concentration of SO2]
  y18 = 0   {caption: SO4; category: Initial concentrations} [Initial concentration of SO4]
  y19 = 0   {caption: NO3; category: Initial concentrations} [Initial concentration of NO3]
  y20 = 0   {caption: N2O5; category: Initial concentrations} [Initial concentration of N2O5]

#output:
   t  {caption: t, min}
  y1  {caption: NO2}
  y2  {caption: NO}
  y3  {caption: O3P}
  y4  {caption: O3}
  y5  {caption: HO2}
  y6  {caption: OH}
  y7  {caption: HCHO}
  y8  {caption: CO}
  y9  {caption: ALD} 
  y10 {caption: MEO2}
  y11 {caption: C2O3}
  y12 {caption: CO2}
  y13 {caption: PAN}
  y14 {caption: CH3O}
  y15 {caption: HNO3}
  y16 {caption: O1D}
  y17 {caption: SO2}
  y18 {caption: SO4}
  y19 {caption: NO3}
  y20 {caption: N2O5}

#parameters:
  k1 = 0.35    {category: Reaction constants} [NO2 -> NO + O3P]
  k2 = 26.6    {category: Reaction constants} [NO + O3 -> NO2]
  k3 = 1.23e4  {category: Reaction constants} [HO2 + NO -> NO2 + OH]
  k4 = 8.6e-4  {category: Reaction constants} [HCHO -> 2 HO2 + CO]
  k5 = 8.2e-4  {category: Reaction constants} [HCHO -> CO]
  k6 = 1.5e4   {category: Reaction constants} [HCHO + OH -> HO2 + CO]
  k7 = 1.3e-4  {category: Reaction constants} [ALD -> MEO2 + HO2 + CO]
  k8 = 2.4e4   {category: Reaction constants} [ALD + OH -> C2O3]
  k9 = 1.65e4  {category: Reaction constants} [C2O3 + NO-> NO2 + MEO2 + CO2]
  k10 = 9e3    {category: Reaction constants} [C2O3 + NO2-> PAN]
  k11 = 0.022  {category: Reaction constants} [PAN-> CH3O + NO2]
  k12 = 1.2e4  {category: Reaction constants} [MEO2 + NO-> CH3O + NO2]
  k13 = 1.88   {category: Reaction constants} [CH3O-> HCHO + HO2]
  k14 = 1.63e4 {category: Reaction constants} [NO2 + OH -> HNO3]
  k15 = 4.8e6  {category: Reaction constants} [O3P -> O3]
  k16 = 3.5e-4 {category: Reaction constants} [O3 -> O1D]
  k17 = 0.0175 {category: Reaction constants} [O3 -> O3P]
  k18 = 1e8    {category: Reaction constants} [O1D -> 2 OH]
  k19 = 4.44e11 {category: Reaction constants} [O1D -> O3P]
  k20 = 1240   {category: Reaction constants} [SO2 + OH -> SO4 + HO2] 
  k21 = 2.1    {category: Reaction constants} [NO3 -> NO]
  k22 = 5.78   {category: Reaction constants} [NO3 -> NO2 + O3P]
  k23 = 0.0474 {category: Reaction constants} [NO2 + O3 -> NO3]
  k24 = 1780   {category: Reaction constants} [NO3 + NO2 -> N2O5]
  k25 = 3.12   {category: Reaction constants} [N2O5 -> NO3 + NO2]

#tolerance: 1e-6`;

const chem = `#name: Chem react
#description: Mass-action kinetics illustration
#comment: 
  Source: https://doi.org/10.1002/ijch.201800003.
#equations:
  dx_1/dt = -k_1 * x_1 + k_2 * (x_2)**2 + k_3 * x_1 * x_3 
           - k_4 * (x_1)**2 - 2 * k_5 * (x_1)**2 + k_6 * x_2 * x_4

  dx_2/dt = 2 * k_1 * x_1 - 2 * k_2 * (x_2)**2 
           + k_5 * (x_1)**2 - k_6 * x_2 * x_4

  dx_3/dt = -k_3 * x_1 * x_3 + k_4 * (x_1)**2 + k_6 * x_2 * x_4

  dx_4/dt = k_5 * (x_1)**2 - k_6 * x_2 * x_4

#parameters:
  k_1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5}
  k_2 = 0.9 {category: Reaction parameters; min: 0.1; max: 5}
  k_3 = 1.2 {category: Reaction parameters; min: 0.1; max: 5}
  k_4 = 3.4 {category: Reaction parameters; min: 0.1; max: 5}
  k_5 = 2.3 {category: Reaction parameters; min: 0.1; max: 5}
  k_6 = 4.5 {category: Reaction parameters; min: 0.1; max: 5}

#inits:
  x_1 = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x_2 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x_3 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x_4 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}

#argument: t
  initial = 0 {units: min; caption: Initial; category: Time; min: 0; max: 5} [Initial time of simulation]
  final = 6 {units: min; caption: Final; category: Time; min: 6; max: 10} [Final time of simulation]
  step = 0.01 {units: min; caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simulation]

#tolerance: 5e-5`;

const adv = `#name: Advanced
#comment:
  This is an advanced template. Modify it. Use multi-line formulas if needed.
  Add new equations, expressions, constants & parameters. Edit these comment lines if required.
#equations:
  dx/dt = E_1 * y + sin(t)

  dy/dt = E_2 * x - pow(pow(t + 1, 5),6)

#expressions:
  E1 = C_1 * exp(-pow(t,2)/(2*C_1)) + P_1
  E2 = C_2 * cos(2 * t) + P_2

#argument: t
  start = 0
  finish = 10 {min: 10; max: 100}
  step = 0.01

#inits:
  x = 2
  y = 0

#constants:
  C_1 = 1
  C_2 = 3

#parameters:
  P_1 = 1
  P_2 = -1

#tolerance: 5e-5`;

const pkpd = `#name: PK-PD
#description: Pharmacokinetic-pharmacodynamic (PK-PD) simulation: two-compartment model
#equations:
  d(depot)/dt = -KA * depot
  d(centr)/dt = KA * depot - CL * C2 - Q * C2 + Q * C3
  d(peri)/dt  = Q * C2 - Q * C3
  d(eff)/dt  = Kin - Kout * (1 - C2/(EC50 + C2)) * eff

#expressions:
  C2 = centr / V2
  C3 = peri / V3

#loop:
  count = 10 {caption: count; category: Dosing; min: 1; max: 20} [Number of doses]
  depot += dose

#argument: t
  start = 0 {units: h; caption: begin; category: Dosing; min: 0; max: 1} [Begin of dosing interval]
  final = 12 {units: h; caption: end; category: Dosing; min: 5; max: 15} [End of dosing interval]
  step = 0.1 {units: h; caption: step; category: Dosing; min: 0.01; max: 0.1} [Time step of simulation]  

#inits:  
  depot = 0 {category: Initial values}
  centr = 0 {category: Initial values} [Central]
  peri = 0 {category: Initial values} [Peripheral]
  eff = 0.2 {category: Initial values} [Effective compartment rate]

#parameters:  
  dose = 1e4 {category: Dosing; min: 1e3; max: 2e4; step: 1e3} [Dosage]
  KA = 0.3 {caption: rate constant; category: Parameters; min: 0.1; max: 1}
  CL = 2 {caption: clearance; category: Parameters; min: 1; max: 5}
  V2 = 4 {caption: central volume; category: Parameters; min: 1; max: 10} [Central compartment volume]
  Q = 1 {caption: inter rate; category: Parameters; min: 0.1; max: 1} [Intercompartmental rate]
  V3 = 30 {caption: peri volume; category: Parameters; min: 20; max: 40} [Peripheral compartment volume]
  EC50 = 8 {caption: effect; category: Parameters; min: 1; max: 10}
  Kin = 0.2 {caption: Kin; category: Parameters; min: 0.1; max: 0.5} [The first-order production constant]
  Kout = 0.2 {caption: Kout; category: Parameters; min: 0.1; max: 0.5} [The first-order dissipation rate constant]
  
#tolerance: 1e-9`;

const rober = `#name: Robertson
#description: Robertson chemical reaction model
#comment: This is classic example of stiff ODEs.
#equations:
  dA/dt = -0.04 * A + 1e4 * B * C
  dB/dt = 0.04 * A - 1e4 * B * C - 3e7 * B**2
  dC/dt = 3e7 * B**2

#inits:
  A = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  B = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  C = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}

#argument: t
  start = 0 {units: sec; caption: Initial; category: Time; min: 0; max: 1} [Initial time of simulation]
  finish = 40 {units: sec; caption: Final; category: Time; min: 2; max: 50} [Final time of simulation]
  step = 0.01 {units: sec; caption: Step; category: Time; min: 0.01; max: 0.1} [Time step of simulation]

#tolerance: 1e-7`;

const ga = `#name: GA-production
#description: Gluconic acid (GA) production by Aspergillus niger modeling
#equations:
  dX/dt = rX
  dS/dt = -gamma * rX - lambda_1 * X
  dO/dt = Kla * (Cod - O) - delta2 * rX - phi * X
  dP/dt = alpha * rX + beta * X

#expressions:
  mu = mu_M * S / (Ks + S) * O / (Ko + O)
  rX = mu * X

#argument: t, 1-st stage
  _t0 = 0 {units: h; caption: initial; category: Misc} [Start of the process]
  _t1 = 60 {units: h; caption: 1-st stage; category: Durations; min: 20; max: 80} [Duration of the 1-st stage]
  step = 0.1 {units: h; caption: step; category: Misc; min: 0.01; max: 1} [Time step of simulation]

#update: 2-nd stage
  duration = overall - _t1
  S += 70

#inits:  
  X = 5 {units: kg/m³; caption: biomass; category: Initial concentrations; min: 1; max: 10} [Aspergillus niger biomass]
  S = 150 {units: kg/m³; caption: glucose; category: Initial concentrations; min: 50; max: 200} [Glucose]
  O = 7 {units: kg/m³; caption: oxygen; category: Initial concentrations; min: 1; max: 10} [Dissolved oxygen]
  P = 0 {units: kg/m³; caption: acid; category: Initial concentrations; min: 0; max: 0.1} [Gluconic acid]

#output:
  t {caption: time}
  X {caption: biomass}
  S {caption: glucose}
  O {caption: oxygen}
  P {caption: acid}

#parameters:
  overall = 100 {units: h; category: Durations; min: 100; max: 140} [Overall duration]
  mu_M = 0.668 {units: 1/h; category: Parameters} [Monod type model parameter]
  alpha = 2.92 {category: Parameters} [Monod type model parameter]
  beta = 0.131 {units: 1/h; category: Parameters} [Monod type model parameter]
  gamma = 2.12 {category: Parameters} [Monod type model parameter]
  lambda = 0.232 {units: 1/h; category: Parameters} [Monod type model parameter]
  delta = 0.278 {category: Parameters} [Monod type model parameter]
  phi = 4.87e-3 {units: 1/h; category: Parameters} [Monod type model parameter]
  Ks = 1.309e2 {units: g/L; category: Parameters} [Monod type model parameter]
  Ko = 3.63e-4 {units: g/L; category: Parameters} [Monod type model parameter]
  Kla = 1.7e-2 {units: 1/s; category: Parameters} [Volumetric mass transfer coefficient]
  Cod = 15 {units: kg/m³; category: Parameters} [Liquid phase dissolved oxygen saturation concentration]
  
#tolerance: 1e-9`;

const model = ga;

const exporter = new ModelToDocExporter(getIVP(model)/*, {eqnTag: FORMULA_TAG.DOUBLE_DOLLAR, valTag: FORMULA_TAG.DOLLAR}*/);
//const exporter = new ModelToLaTeXExporter(getIVP(model));

//console.log(exporter.getDocLines().join('\n\n'));
console.log(exporter.getShortDocLines().join('\n\n'));
