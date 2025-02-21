import * as DGL from "./index"; 

const model = `#name: UDel
#tags: model
#description: Bioreactor simulation, ver. Nov 19, 2024

#equations:
  d(X_T)/dt  =   mu * X_V - K_lysis * (X_T - X_V)
  d(X_V)/dt  =  (mu - k_d) * X_V
  d(Glc)/dt  =   dU3
  d(Gln)/dt  =   dU4
  d(Lac)/dt  =   dU5
  d(Amm)/dt  =  - Y_Amm_Gln * dU4
  d(ATP)/dt  =   dU7
  d(MAb)/dt  =   dMab
  d(Vol)/dt  =   0

#expressions:
  mu_Glc        = Glc / (K_Glc + Glc)
  mu_Gln        = Gln / (K_Gln + Gln)
  Glc_per_Cell  = Glc / X_V
  sigma_Glc     = 1 / (1 + exp(sigm_factor * (Glc_per_Cell - Glc_Cell_Th)))
  mu_Lac        = mu_1_max * sigma_Glc * K_Lac_i / (K_Lac_i + Lac) + mu_2_max * (1 - sigma_Glc) * Lac / (K_Lac + Lac)
  mu_Amm        = K_Amm / (K_Amm + Amm)
  mu            = mu_Gln * mu_Amm * mu_Lac * mu_Glc
  k_d           = k_d_max * k_mu / (mu + k_mu)
  dU3           = - mu / Y_X_Glc * X_V - m_Glc * X_V
  dU4           = - mu / Y_X_Gln * X_V - m_Gln * X_V
  dU5           = - Y_Lac_Glc * dU3 * sigma_Glc - mu / Y_X_Lac * X_V * (1 - sigma_Glc)
  dU7           = - sigma_Glc * 30 * dU3 - (1 - sigma_Glc) * 12 * dU5 - 15 * dU4
  dMab          = (1/1342 * dU7 - 1/66 * dU3 - 1/62 * dU4)
  Prod          = dMab
  denom         = X_V * W_dry_scale * W_dry
  SPr           = dMab / denom
  SSUR          = -dU3 / denom

#argument: t
  initial = 0     {units: h; category: Time; min: 0; max: 10}   [Begin of simulation interval]
  final = 250   {units: h; category: Time; min: 20; max: 750} [Feeding interval interval]
  step = 3      {units: h; category: Time; min: 0.5; max: 24} [Time step of simulation]

#output:
  t
  X_T
  X_V
  Glc
  Gln
  Lac
  Amm
  MAb
  ATP
  Vol
  Prod
  SPr
  SSUR

#inits:
  X_T = 1      {category: Initials; units: mln. cells / L}   [Initial total cell density]
  X_V = 1      {category: Initials; units: mln. cells / L}   [Initial viable cell density]
  Glc = 1      {category: Initials; units: g/L}              [Initial glucose concentration]
  Gln = 1      {category: Initials; units: g/L}              [Initial glutamine concentration]
  Lac = 0.5    {category: Initials; units: g/L}              [Initial lactate concentration]
  Amm = 1.1    {category: Initials; units: g/L}              [Initial ammonium concentration]
  MAb = 0      {category: Initials; units: g/L}              [Initial antibody concentration]
  ATP = 0      {category: Initials; units: g/L}              [Initial ATP concentration]
  Vol = 2      {category: Initials; units: L; format: #0.00} [Initial volume] 

#parameters:
  mu_1_max        = 0.048     {caption:  mu_1_max ; category:  Parameters ; units:  1/h ; min: 0.00000001; max: 1500; step: 0.00001} [Maximum growth rate (exponential)]
  mu_2_max        = 0.012     {caption:  mu_2_max ; category:  Parameters ; units:  1/h ; min: 0.00000001; max: 1500; step: 0.00001} [Maximum growth rate (stationary)]
  K_Glc           = 1         {caption:  K_Glc ; category:  Parameters ; units:  g/L ; min: 0.00000001; max: 10; step: 0.0000001} [Monod constant of glucose]
  K_Gln           = 0.22      {caption:  K_Gln ; category:  Parameters ; units:  g/L ; min: 0.00000001; max: 10; step: 0.0000001} [Monod constant of glutamine]
  K_Lac           = 0.2       {caption:  K_Lac ; category:  Parameters ; units:  g/L ; min: 0.00000001; max: 0.3; step: 0.0000001} [Monod constant of lactate]
  K_Lac_i         = 150       {caption:  K_Lac_i ; category:  Parameters ; units: g/L ; min: 0.00000001; max: 30000.0; step: 1.0} [Constant of lactate inhibition]
  K_Amm           = 40        {category:  Parameters}
  k_mu            = 0.01      {category:  Parameters}
  k_d_max         = 0.01      {category:  Parameters}
  m_Glc           = 8e-6      {category:  Parameters; format: #0.0E00}
  m_Gln           = 0         {category:  Parameters}
  Y_X_Glc         = 15        {category:  Parameters}
  Y_X_Gln         = 190       {category:  Parameters}
  Y_Lac_Glc       = 4.0       {category:  Parameters}
  Y_X_Lac         = 7.500075  {category:  Parameters; min: 0.0001; max: 30}
  Y_Amm_Gln       = 1.1       {category:  Parameters}
  K_lysis         = 0.000001  {category:  Parameters}
  Glc_Cell_Th     = 0.2       {category:  Parameters; min: 0.001; max: 0.4}
  sigm_factor     = -10       {category:  Parameters; min: -100; max: -10}
  W_dry           = 2.5e-12   {category:  Constants; units: g/cell; format: #0.0E00}

#constants:
  W_dry_scale = 1e7

#tolerance: 0.0000001`;

const ivp = DGL.getIVP(model);

//console.log(ivp);

const ivpWW = DGL.getIvp2WebWorker(ivp);

console.log(ivpWW);

const inputs = new Float64Array([
    0,
    250,
    10,
    1,
    1,
    1,
    1,
    0.5,
    1.1,
    0,
    0,
    2,
    0.048,
    0.012,
    1,
    0.22,
    0.2,
    150,
    40,
    0.01,
    0.01,
    0.000008,
    0,
    15,
    190,
    4,
    7.500075,
    1.1,
    0.000001,
    0.2,
    -10,
    2.50E-12,
]);

const pipeline: DGL.Pipeline = {
    wrappers: [
        {
            preproc: null,
            postproc: null,
        },
        {
            preproc: null,
            postproc: null,
        },
    ],
    out: null,
};

const solution = DGL.applyPipeline(pipeline, ivpWW, inputs);

const length = solution[0].length;

for (let i = 0; i < length; ++i)
    console.log(solution[0][i], '    ', solution[3][i], '  ', solution[4][i]);

console.log(solution);
