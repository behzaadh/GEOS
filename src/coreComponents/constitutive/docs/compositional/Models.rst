List of available compositional fluid models
=============================================

The following table outlines the available compositional fluid models (specified via their XML elements), detailing the number of phases they support, their underlying thermodynamic phase split (flash) methods, the phase density and viscosity models used, and their recommended applications:

.. list-table:: 
   :widths: 100
   :header-rows: 1

   * - Compositional fluid models
   * - ``CompositionalTwoPhaseFluid`` :ref:`XML_CompositionalTwoPhaseFluid`

       * Number of phases: 2.
       * Possible applications: Idealised simulations, synthetic testing, or when fluid flow is dominated entirely by pressure gradients.
       * Phase split method: Negative Two-Phase Flash.
       * Phase density models: Compositional density.
       * Phase viscosity models: Constant viscosity.

   * - ``CompositionalTwoPhaseFluidLohrenzBrayClark`` :ref:`XML_CompositionalTwoPhaseFluidLohrenzBrayClark`

       * Number of phases: 2.
       * Possible applications: Miscible gas injection, primary depletion of volatile oils, and scenarios where fluid compositions change drastically.
       * Phase split method: Negative Two-Phase Flash.
       * Phase density models: Compositional density.
       * Phase viscosity models: Lohrenz-Bray-Clark viscosity.

   * - ``CompositionalTwoPhaseFluidPhillipsBrine`` :ref:`XML_CompositionalTwoPhaseFluidPhillipsBrine`

       * Number of phases: 2.
       * Possible applications: Carbon sequestration (CO2-brine systems) or water-flooding in oil reservoirs.
       * Phase split method: Negative Two-Phase Flash.
       * Phase density models: Phillips brine density model for the aqueous phase; Compositional density for the gas phase.
       * Phase viscosity models: Phillips brine viscosity model for the aqueous phase; Lohrenz-Bray-Clark viscosity for the gas phase.

   * - ``CompositionalThreePhaseFluidLohrenzBrayClark`` :ref:`XML_CompositionalThreePhaseFluidLohrenzBrayClark`

       * Number of phases: 3.
       * Possible applications: Multi-phase reservoir simulations where water solubility in the hydrocarbon phases is negligible.
       * Phase split method: Immiscible Water Flash.
       * Phase density models: Immiscible water density for the aqueous phase; Compositional density for the hydrocarbon phases.
       * Phase viscosity models: Immiscible water viscosity for the aqueous phase; Lohrenz-Bray-Clark viscosity for the hydrocarbon phases.

   * - ``CompositionalTwoPhaseKValueFluidLohrenzBrayClark`` :ref:`XML_CompositionalTwoPhaseKValueFluidLohrenzBrayClark`

       * Number of phases: 2.
       * Possible applications: Thermal recovery processes or relatively stable black-oil systems.
       * Phase split method: K-Value Flash.
       * Phase density models: Compositional density.
       * Phase viscosity models: Lohrenz-Bray-Clark viscosity.

   * - ``CompositionalTwoPhaseKValueFluidPhillipsBrine`` :ref:`XML_CompositionalTwoPhaseKValueFluidPhillipsBrine`

       * Number of phases: 2.
       * Possible applications: Saline aquifers and CO2 sequestration with stable compositions where pre-tabulated K-values are sufficient.
       * Phase split method: K-Value Flash.
       * Phase density models: Phillips brine density model; Compositional density.
       * Phase viscosity models: Phillips brine viscosity model; Lohrenz-Bray-Clark viscosity.

   * - ``CompositionalTwoPhaseKValueThermalFluidLohrenzBrayClark`` :ref:`XML_CompositionalTwoPhaseKValueThermalFluidLohrenzBrayClark`

       * Number of phases: 2.
       * Possible applications: Thermal oil–gas recovery with tabulated K(P,T) and compositional enthalpy.
       * Phase split method: Two-phase K-Value Flash.
       * Phase density models: Compositional density.
       * Phase viscosity models: Lohrenz-Bray-Clark viscosity.
       * Phase enthalpy model: Compositional enthalpy (Poling heat capacity plus cubic departure).

   * - ``CompositionalTwoPhaseKValueThermalFluidPhillipsBrine`` :ref:`XML_CompositionalTwoPhaseKValueThermalFluidPhillipsBrine`

       * Number of phases: 2.
       * Possible applications: Thermal gas–water systems with tabulated K(P,T), including CO2 dissolution in brine.
       * Phase split method: Two-phase K-Value Flash.
       * Phase density models: Phillips brine density for the aqueous phase; Compositional density for the gas phase.
       * Phase viscosity models: Phillips brine viscosity for the aqueous phase; Lohrenz-Bray-Clark viscosity for the gas phase.
       * Phase enthalpy model: Compositional enthalpy.

   * - ``CompositionalThreePhaseKValueThermalFluidLohrenzBrayClark`` :ref:`XML_CompositionalThreePhaseKValueThermalFluidLohrenzBrayClark`

       * Number of phases: 3.
       * Possible applications: Thermal oil/gas/water flow when CO2 (or other components) must dissolve in water according to tabulated aqueous K-values.
       * Phase split method: Three-phase K-Value Flash (vapour–liquid and aqueous–liquid K-tables).
       * Phase density models: Compositional density for oil, gas, and aqueous phases.
       * Phase viscosity models: Lohrenz-Bray-Clark viscosity for all phases.
       * Phase enthalpy model: Compositional enthalpy.