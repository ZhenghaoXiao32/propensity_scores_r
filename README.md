# Propensity score analysis using R


## Introduction
### Rubin's Causal Model

We start with  Rubin's casual model. The causal effect of a treatment (<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Z_{i}=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Z_{i}=1" title="Z_{i}=1" /></a>) over the control (<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Z_{i}=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Z_{i}=0" title="Z_{i}=0" /></a>), for a given individual or unit (<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y_{i}" title="Y_{i}" /></a>) and an interval of time from 0 to 1 (<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{0}_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{0}_{i}" title="Y^{0}_{i}" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}" title="Y^{1}_{i}" /></a>) is the ***difference** between **what would have happened at time 1 if the unit had been exposed to treatment initiated at time 0 (<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}|Z_{i}=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}|Z_{i}=1" title="Y^{1}_{i}|Z_{i}=1" /></a>)** and **what would have happened at time 1 if the unit had been exposed to control initiated at time 0 (<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}|Z_{i}=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}|Z_{i}=0" title="Y^{1}_{i}|Z_{i}=0" /></a>)**.*

In Rubin's causal model, we have in fact four potential outcomes: 
* at time 0 with control:<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{0}_{i}|Z_{i}=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{0}_{i}|Z_{i}=0" title="Y^{0}_{i}|Z_{i}=0" /></a>
* at time 0 with treatment:<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{0}_{i}|Z_{i}=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{0}_{i}|Z_{i}=1" title="Y^{0}_{i}|Z_{i}=1" /></a>
* at time 1 with control:<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}|Z_{i}=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}|Z_{i}=0" title="Y^{1}_{i}|Z_{i}=0" /></a>
* at time 1 with treatment:<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}|Z_{i}=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}|Z_{i}=1" title="Y^{1}_{i}|Z_{i}=1" /></a>

For example, we have Tom to test if a new drug can reduce his blood pressure level or not. We tested his blood pressure at time 0, so we got <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{0}_{i}|Z_{i}=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{0}_{i}|Z_{i}=0" title="Y^{0}_{i}|Z_{i}=0" /></a>, then we asked Tom to take the pill, and tested his blood pressure again after an hour, now we got <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}|Z_{i}=1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}|Z_{i}=1" title="Y^{1}_{i}|Z_{i}=1" /></a>. But we what have now is not enough to figure out the causal effect, we need <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y^{1}_{i}|Z_{i}=0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;Y^{1}_{i}|Z_{i}=0" title="Y^{1}_{i}|Z_{i}=0" /></a> which is the blood pressure level of Tom at time 1 if he had not taken the pill (or had taken some kind of placebo).

You should have noticed that in the example what would have happend at time 1 if the unit had been exposed to control initiated at time 0 is impossible to observe for a individual who had been exposed to treatment initiated at time 0. It is obvious that we cannot know the both results of two different conditions of the same time. However, they are what we needed to calculate the causal effect of a treatment. That is called the fundamental problem of causal inference. 

### Stable Unit Treatment Value Assumption

We need an assumption to solve the fundamental problem. The assumption is: "The potential outcome of the observation on one unit should be unaffected by the particular assignment of treatments to the other unit" (Cox 1958). The math formula for it is: given observed covariates <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;X:\left(&space;Y^{0},Y^{1}\right)&space;\bot&space;Z|X" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;X:\left(&space;Y^{0},Y^{1}\right)&space;\bot&space;Z|X" title="X:\left( Y^{0},Y^{1}\right) \bot Z|X" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;0&space;<&space;p\left(&space;Z_{i}^{1}|X\right)&space;<&space;1" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;0&space;<&space;p\left(&space;Z_{i}^{1}|X\right)&space;<&space;1" title="0 < p\left( Z_{i}^{1}|X\right) < 1" /></a> (Rosenbaum & Rubin, 1983).

We can use a example to better illustrate it: Tom's change in blood pressure may depend on whether or not Jerry receives the drug. Assuming Tom and Jerry live together and Jerry always cooks. The drug causes change of Jerry's blood pressure and makes Jerry want to eat salty foods. So Jerry cooks food with more salt and it increases Tom's blood pressure. In that situation, we will say Tom's potential outcome is correlated to which treatment Jerry receives. That is a volation of the stable unit treatment value assumption. 

### How We Solve The Problem with Experimental Design

To solve the fundamental problem of causal inference, we have two methods in essence. The first one is to repeat the experiment with the same individual or unit (repeated measure). As it is the same individual or unit, the main noise should be associated with different time. Generally we need to assume that responses of the individual, we say Tom, is independent of each other. This assumption might be unstable in real life, like if the blood concentratioin of the drug needs some time (maybe a month) to decrease, and the interval of the measurements is less than that time, then we violates the assumption. However, in experimental design, we can a particular design to cancel the bias out:

month | first month | second month | third month | fourth month
----- | ----------- | ------------ | ----------- | ------------
last month | control | treatment | treatment | control
this month | treatment | treatment | control | control

The problem of this design is it increases the number of treatment assignment and makes the calculation of average causal effect harder. Another method is multiple subjects (by groups). We can recruit multiple people for the hypertension drug test, assign them to treatment group and control group. The difference of different individuals can be reduced by random assignment.
However, we are not able to control everything, in most studies, we cannot use random assignment nor manipulation of conditions. Research designs to estimate treatment effects but have no random assignment to conditions is called quasi-experimental or nonexperimental designs. The propensity scores and methods of apply them are designed to reduce the bias the quasi-experimental designs.

### Types of Treatment Effects

Before we get started with propensity scores, we still need to know about the types of treatment effects:
* The average treatment effect (ATE) is the difference between the expected values of the potential outcomes of all individuals in the treated and untreated conditions: <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;ATE=E\left(&space;Y_{i}^{1}\right)&space;-E\left(&space;Y_{i}^{0}\right)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;ATE=E\left(&space;Y_{i}^{1}\right)&space;-E\left(&space;Y_{i}^{0}\right)" title="ATE=E\left( Y_{i}^{1}\right) -E\left( Y_{i}^{0}\right)" /></a>
* The average treatment effect on the treated (ATT) is the difference between the expected values of the potential outcomes of treated individuals: <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;ATT=E\left(&space;Y_{i}^{1}|Z=1\right)&space;-E\left(&space;Y_{i}^{0}|Z=1\right)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;ATT=E\left(&space;Y_{i}^{1}|Z=1\right)&space;-E\left(&space;Y_{i}^{0}|Z=1\right)" title="ATT=E\left( Y_{i}^{1}|Z=1\right) -E\left( Y_{i}^{0}|Z=1\right)" /></a>
* The average treatment effect on the untreated (control) (ATC) is the difference between the expected values of the potential outcomes of the untreated individuals: <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;ATC=E\left(&space;Y_{i}^{1}|Z=0\right)&space;-E\left(&space;Y_{i}^{0}|Z=0\right)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\inline&space;ATC=E\left(&space;Y_{i}^{1}|Z=0\right)&space;-E\left(&space;Y_{i}^{0}|Z=0\right)" title="ATC=E\left( Y_{i}^{1}|Z=0\right) -E\left( Y_{i}^{0}|Z=0\right)" /></a>

In experimental designs, the ATE is euqal to the ATT and the ATC because of the balance achieved by random assignment. But in quasi-experimental design, these values could differ substantially. So, what type of treatment effect we can research in quasi-experimental design is depended on whether we achieved balance in the particular group.


## Contents
- [x] Propensity scores weighting for ATT
- [x] Propensity scores weighting for ATE
- [x] Propensity scores stratification for ATT
- [ ] Propensity scores stratification for ATE
- [ ] Propensity scores matching for ATT
- [ ] Propensity scores matching for ATE

## Maintenances
- [x] Fix unclear variable names in stratification part
