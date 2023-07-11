## Overview

<p align="center">
  <img width="980" height="600" src="https://github.com/UM-Myeloma-Genomics/GCP_MM/blob/main/guide/figs/git_overview.png?raw=true">
</p>

## User guide for website

The website can be accessed at https://irmma-risk-calculator.miami.edu

For a detailed breakdown of the dataset across various states included in our knowledge bank, please refer to Figure 3a in the paper.

For query suggestions based on what IRMMa has seen and learned from,
<p align="center">
  <img width="600" height="500" src="./figs/recommendations.png">
</p>

### IRMMa Risk Calculator manual

This website is a tool to calculate risks according to an individual’s characteristics. 

* Please first complete the Attribute Form which is the set of features used to calculate multistate/treatment risks. <br /> *Age* (*Demographics*), *ISS*, and Induction (Phase 1) (*Therapy*) are **mandatory fields**. <br /> If there is no selection on *IGH-Translocations* or *Genomics*, the default is NA. 
* When you have completed your selections, press the *Submit* button. Of note, the threshold for  minimum number of cases that IRMMa utilizes to make a prediction is 1.

**Once submitted, IRMMa calculates risks and the following is displayed:**

* **State Risk** - Multistate risks from start of P1 up until 5 years.
* **Treatment-Adjusted Risk** - Risk of POD and/or death at 5 years based on selected and alternative P1 treatment combinations and P2 treatments, MAINT/CONT.TREAT and HDM-ASCT.

**For interpretation of risks**, please refer to figure below, where the risks are of the query with respect to the knowledge bank scores.
<p align="center">
  <img width="720" height="300" src="./figs/risks.png">
</p>

## Example usage

![alt-text-1](./figs/ex1_2.png "High risk") ![alt-text-2](./figs/ex1_3.png "Low risk")


**Acronyms/definitions**:

* **IRMMa** - Individualized Risk Model for Myeloma.
* **MAINT/CONT.TREAT (M.C.T.)** - Maintenance/continuous treatment (Defined as either maintenance therapy continuing after induction or a treatment lasting more than 12 months, Palumbo et al. JCO 2015).
* **CHEMO** - Any chemotherapy regiment received except high-dose melphalan. This includes cyclophosphamide, low-dose melphalan, or any platinum-based chemotherapy.
* **HDM-ASCT** - High-dose melphalan followed by autologous stem cell transplantation.
* **Induction, Phase 1 (P1)** - Induction treatment within 1 year.
* **Post-Induction, Phase 2 (P2)** - Addition of HDM-ASCT and/or MAINT/CONT.TREAT or neither upon completion of P1 (without POD).
* **POD** - Progression of Disease.
* **Probability of being alive** - Probability of *Alive in P1 + Alive after POD (P1) + Alive in P2 + Alive after POD (P2)*.
* **Probability of POD** - Probability of *Alive after POD (P1) + Alive after POD (P2) + Death after POD (P1) + Death after POD (P2)*. 
* **Risk of POD and/or deceased** - Probability of *Alive in P1 + Alive after POD (P1) + Death after POD (P1) + Death in P2 + Alive after POD (P2) + Death after POD (P2)*.
* In **Treatment-Adjusted Risk** plot, for each therapy $i$ where $i$ are all the combination therapies that is present in our knowledge bank, we provide the variance **var (%)** as the difference in risk between query ${q}$ selected therapy in P1, P2 and alternative therapy $i$. $var = {risk_i - risk_q \over risk_q}\times100$. Thus negative var indicates reduced risk of POD and/or death in comparison to query selected therapy, positive indicates worse risk of POD and/or death, 0 indicates risk of query $q$. 

**More Information**
* [Data preparation and genomic classification](https://github.com/UM-Myeloma-Genomics/GCP_MM/tree/main/genomic)
* [IRMMa](https://github.com/UM-Myeloma-Genomics/GCP_MM/tree/main/prognostication)

