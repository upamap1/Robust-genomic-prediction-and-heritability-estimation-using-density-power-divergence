{\rtf1\ansi\ansicpg1252\cocoartf2709
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Install required packages for robust genomic prediction\
\
required_packages <- c(\
  "lme4",\
  "nlme",\
  "robustlmm",\
  "rrBLUP",\
  "Metrics",\
  "psych",\
  "matrixStats",\
  "insight",\
  "caTools",\
  "caret",\
  "glmnet",\
  "ie2misc",\
  "foreach",\
  "doParallel",\
  "MASS",\
  "Agricolae",\
  "AlphaSimR"\
)\
\
installed <- rownames(installed.packages())\
for (pkg in required_packages) \{\
  if (!pkg %in% installed) \{\
    install.packages(pkg, dependencies = TRUE)\
  \}\
\}}