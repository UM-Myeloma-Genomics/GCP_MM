## User guide for website

The website can be accessed [here](https://irmma-risk-calculator.miami.edu).

For dataset breakdown in the different states present in our knowledge bank is shown in  Fig 3a in the paper.

Please refer to [recommendations](recommendations.pdf) on what queries might provide back results by IRMMa.

### IRMMa Risk Calculator manual

This website is a tool to calculate risks according to an individual’s characteristics. The Attribute Form is the set of features used to calculate multistate/treatment risks. *Age* (*Demographics*), *ISS*, and Induction (Phase 1) (*Therapy*) are mandatory fields. Please select them in the section to the left by clicking one of the section drop-downs. If there is no selection on *IGH-Translocations* or *Genomics*, the default is NA. When you have completed your selections, press the Submit button.

Once submitted, IRMMa calculates risks and the following is displayed:

* **State Risk** - Multistate risks from start of P1 up until 5 years.
* **Treatment-Adjusted Risk** - Risk of POD or death at 5 years based on selected and alternative first-line treatment combinations and P2 treatments, MAINT/CONT.TREAT and HDM-ASCT.

**Acronyms/definitions**:

* **IRMMa** - Individualized Risk Model for Myeloma.
* **MAINT/CONT.TREAT (M.C.T.)** - Maintenance/continuous treatment (Defined as either maintenance therapy continuing after induction or a treatment lasting more than 12 months, Palumbo et al. JCO 2015).
* **CHEMO** - Any chemotherapy regiment received except high-dose melphalan. This includes cyclophosphamide, low-dose melphalan, or any platinum-based chemotherapy.
* **HDM-ASCT** - High-dose melphalan followed by autologous stem cell transplantation.
* **Induction, Phase 1 (P1)** - Induction treatment within 1 year.
* **Consolidation, Phase 2 (P2)** - Addition of HDM-ASCT and/or MAINT/CONT.TREAT or neither upon completion of P1 (without POD).
* **POD** - Progression of Disease.

**More Information**
* [Data preparation and genomic classification](https://github.com/UM-Myeloma-Genomics/GCP_MM/tree/main/genomic)
* [IRMMa](https://github.com/UM-Myeloma-Genomics/GCP_MM/tree/main/prognostication)
