{\rtf1\ansi\ansicpg936\cocoartf2822
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Figure 1 PyMOL colouring script\
# Object names should match the loaded PDB object names:\
# fold_cd276_monomer_model_0, antibody_input, nanobody_input,\
# IgG_FP_0037, NB_FP_0087\
\
hide everything\
show cartoon\
bg_color white\
\
# Improve visual clarity\
set cartoon_fancy_helices, 1\
set antialias, 2\
set ray_opaque_background, off\
\
# -----------------------------\
# Panel A: CD276 monomer\
# -----------------------------\
color green, fold_cd276_monomer_model_0 and chain A\
color red, fold_cd276_monomer_model_0 and chain A and resi 242-248\
color yellow, fold_cd276_monomer_model_0 and chain A and resi 282-285\
\
# -----------------------------\
# Panel B/D: IgG antibody complex\
# antibody_input = AF3-repredicted complex or input complex, depending on panel\
# IgG_FP_0037 = PPIFlow/Flowpacker output\
# -----------------------------\
\
# AF3 antibody input / re-predicted structure\
color green, antibody_input and chain A\
color red, antibody_input and chain A and resi 242-248\
color yellow, antibody_input and chain A and resi 282-285\
color marine, antibody_input and chain B\
color cyan, antibody_input and chain C\
\
# PPIFlow IgG output\
color green, IgG_FP_0037 and chain C\
color red, IgG_FP_0037 and chain C and resi 242-248\
color yellow, IgG_FP_0037 and chain C and resi 282-285\
color marine, IgG_FP_0037 and chain A\
color cyan, IgG_FP_0037 and chain B\
\
# -----------------------------\
# Panel C/E: Nanobody complex\
# nanobody_input = AF3-repredicted complex or input complex, depending on panel\
# NB_FP_0087 = PPIFlow/Flowpacker output\
# -----------------------------\
\
# AF3 nanobody input / re-predicted structure\
color green, nanobody_input and chain A\
color red, nanobody_input and chain A and resi 242-248\
color yellow, nanobody_input and chain A and resi 282-285\
color magenta, nanobody_input and chain B\
\
# PPIFlow nanobody output\
color green, NB_FP_0087 and chain C\
color red, NB_FP_0087 and chain C and resi 242-248\
color yellow, NB_FP_0087 and chain C and resi 282-285\
color magenta, NB_FP_0087 and chain A\
\
# Highlight selected epitope residues\
show sticks, (resi 242-248 or resi 282-285)\
set stick_radius, 0.18}