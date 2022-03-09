
# This script calculates a pH given a lifetime (assuming a 4 parameter
# logistic fit model for mScarlet lifetime)

def pH_from_tau(tau, pKa, min_val, max_val, hill):
    return pKa * ((min_val - max_val)/(tau-max_val) - 1)**(1/hill)

def tau_from_pH(pH, pKa, min_val, max_val, hill):
    return max_val + (min_val - max_val) / (1 + (pH/pKa)**hill)

# U2OS lysosome calibration
pKa = 4.90
min_val = 1.72
max_val = 3.25
hill = 13

tau = 2.71

print(round(pH_from_tau(tau, pKa, min_val, max_val, hill), 2))

pH = 6.3

print(round(tau_from_pH(pH, pKa, min_val, max_val, hill), 2))
