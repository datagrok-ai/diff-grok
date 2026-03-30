/** Example: pollution model (20 species) → LaTeX, no metadata, no inits. */

import {convertIvpToLatex} from '../index';

const MODEL = `#name: Pollution
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
  t0 = 0   {units: min; caption: Initial; category: Time; min: 0; max: 0.9}
  t1 = 60  {units: min; caption: Final; category: Time; min: 1; max: 100; step: 1}
  h = 0.1  {units: min; caption: Step; category: Time; min: 0.001; max: 0.1; step: 0.001}

#inits:
  y1 = 0    {caption: NO2; category: Initial concentrations; min: 0; max: 0.1}
  y2 = 0.2  {caption: NO; category: Initial concentrations; min: 0; max: 0.4}
  y3 = 0    {caption: O3P; category: Initial concentrations; min: 0; max: 0.1}
  y4 = 0.04 {caption: O3; category: Initial concentrations; min: 0; max: 0.1}
  y5 = 0    {caption: HO2; category: Initial concentrations}
  y6 = 0    {caption: OH; category: Initial concentrations}
  y7 = 0.1  {caption: HCHO; category: Initial concentrations}
  y8 = 0.3  {caption: CO; category: Initial concentrations}
  y9 = 0.01 {caption: ALD; category: Initial concentrations}
  y10 = 0   {caption: MEO2; category: Initial concentrations}
  y11 = 0   {caption: C2O3; category: Initial concentrations}
  y12 = 0   {caption: CO2; category: Initial concentrations}
  y13 = 0   {caption: PAN; category: Initial concentrations}
  y14 = 0   {caption: CH3O; category: Initial concentrations}
  y15 = 0   {caption: HNO3; category: Initial concentrations}
  y16 = 0   {caption: O1D; category: Initial concentrations}
  y17 = 0.007 {caption: SO2; category: Initial concentrations}
  y18 = 0   {caption: SO4; category: Initial concentrations}
  y19 = 0   {caption: NO3; category: Initial concentrations}
  y20 = 0   {caption: N2O5; category: Initial concentrations}

#output:
   t  {caption: t, min}
  y1  {caption: NO2}
  y2  {caption: NO}

#parameters:
  k1 = 0.35    {category: Reaction constants}
  k2 = 26.6    {category: Reaction constants}
  k3 = 1.23e4  {category: Reaction constants}
  k4 = 8.6e-4  {category: Reaction constants}
  k5 = 8.2e-4  {category: Reaction constants}
  k6 = 1.5e4   {category: Reaction constants}
  k7 = 1.3e-4  {category: Reaction constants}
  k8 = 2.4e4   {category: Reaction constants}
  k9 = 1.65e4  {category: Reaction constants}
  k10 = 9e3    {category: Reaction constants}
  k11 = 0.022  {category: Reaction constants}
  k12 = 1.2e4  {category: Reaction constants}
  k13 = 1.88   {category: Reaction constants}
  k14 = 1.63e4 {category: Reaction constants}
  k15 = 4.8e6  {category: Reaction constants}
  k16 = 3.5e-4 {category: Reaction constants}
  k17 = 0.0175 {category: Reaction constants}
  k18 = 1e8    {category: Reaction constants}
  k19 = 4.44e11 {category: Reaction constants}
  k20 = 1240   {category: Reaction constants}
  k21 = 2.1    {category: Reaction constants}
  k22 = 5.78   {category: Reaction constants}
  k23 = 0.0474 {category: Reaction constants}
  k24 = 1780   {category: Reaction constants}
  k25 = 3.12   {category: Reaction constants}

#tolerance: 1e-6`;

// Large model: no metadata, no inits — focus on equations & parameters
const latex = convertIvpToLatex(MODEL, {
  includeMetadata: false,
  includeInits: false,
});
console.log(latex);
