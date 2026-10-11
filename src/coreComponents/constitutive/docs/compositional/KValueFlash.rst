K-value flash
-------------

The K-value flash model offers a computationally simplified alternative to the rigorous equation of state flash calculation. Rather than iteratively evaluating complex equations of state to find fugacity coefficients, this model determines phase splits using pre-tabulated K-values. 

The model supports both two-phase and three-phase variants. The two-phase flash calculates equilibrium partitioning between a liquid and a vapour phase. The three-phase flash introduces an additional aqueous phase. Unlike the immiscible-water flash, the three-phase K-value flash **allows non-water components (for example CO2) to dissolve in the aqueous phase**. The dissolved mole fractions are whatever the user tabulates in the aqueous K-values :math:`K_i^W(P,T)`.

Methodology
~~~~~~~~~~~

In this formulation, the equilibrium ratio for each component in each secondary phase relative to the primary liquid phase is treated strictly as a predefined function of pressure and temperature:

* Vapour vs liquid: :math:`K_i^V = y_i / x_i`
* Aqueous vs liquid (three-phase only): :math:`K_i^W = w_i / x_i`

Users input these functions as either tables or symbolic functions. To optimize performance during the simulation, the model discretizes these functions during initialization to construct a multi-dimensional hypercube of K-values spanning the phase, component, pressure, and temperature dimensions.

The pressure and temperature coordinates for this hypercube grid can be explicitly chosen and specified by the user. If they are omitted, the model automatically inspects the provided table functions, extracts all unique coordinate points from the underlying data, and constructs the grid dynamically. If only a single point is found, an artificial offset is added to permit interpolation.

``kValueTables`` order
~~~~~~~~~~~~~~~~~~~~~~

For an :math:`N`-phase, :math:`N_c`-component system the list has :math:`(N-1)\times N_c` function names, stored with **secondary phase index outermost and component index innermost** (all vapour–liquid K-tables, then all aqueous–liquid K-tables):

.. math::

   \{ K_1^V,\ldots,K_{N_c}^V,\; K_1^W,\ldots,K_{N_c}^W \}

For a four-component oil/gas/water mixture ``{CO2, C1C2, C4+, H2O}`` that is eight tables. Water must be a listed component (dead-water / immiscible-water models cannot represent CO2 in brine). Hydrocarbons that should stay out of water are given small :math:`K_i^W`; CO2 dissolution is controlled by :math:`K_{\mathrm{CO2}}^W(P,T)`. Those aqueous tables are generated **offline** (Intersect, CoolProp, experiments) the same way oil–gas K-tables are built.

Mapping Intersect-style oil–gas K plus a user aqueous table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intersect decks such as ``T_IXD14`` typically tabulate only oil–gas K vs pressure at several temperatures, and treat water as a dead phase. To run the same system in GEOS with CO2–water solubility:

* Convert each Intersect oil–gas ``KValuePressureTable`` slice to a GEOS 2D ``TableFunction`` of ``(pressure, temperature)`` for every hydrocarbon (and for H2O vapour–liquid K).
* Supply **new** aqueous–liquid tables: a large :math:`K_{\mathrm{H2O}}^W`, a solubility-based :math:`K_{\mathrm{CO2}}^W`, and small values for C1C2 / C4+.
* For thermal runs add ``componentHeatCapacityCoefficients`` (Poling 5-coefficient Cp) and optional ``referenceEnthalpy`` / ``enthalpyReferenceTemperature``. Use a thermal catalog name such as ``CompositionalThreePhaseKValueThermalFluidLohrenzBrayClark`` and set ``isThermal="1"`` on the flow solver and ``WellManager``.

Equilibrium calculation
~~~~~~~~~~~~~~~~~~~~~~~

During the simulation step, the model performs bilinear interpolation on the pre-built hypercube to evaluate the K-values and their partial derivatives at the current pressure and temperature.

Because the K-values are explicitly known, the iterative fugacity updates are bypassed. For two phases the vapour fraction :math:`V` is found by solving the Rachford-Rice equation (Rachford and Rice, 1952) using successive substitution and Newton iterations:

.. math::
    \sum_{i=1}^{N_c} \frac{z_i (K_i - 1)}{1 + V(K_i - 1)} = 0

Once the vapour fraction is determined, the phase compositions are directly evaluated:

.. math::
    x_i = \frac{z_i}{1 + V(K_i - 1)}, \quad y_i = K_i x_i

This direct solution approach intrinsically supports a "negative" flash. The mathematical solver allows the vapour fraction, :math:`V`, to temporarily converge to values outside the physically meaningful bounds of :math:`[0, 1]`. If the resulting :math:`V \le 0` or :math:`V \ge 1`, the system cleanly truncates the compositions to represent a single-phase liquid or single-phase vapour solution, respectively.

Three-phase Rachford–Rice
~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`L = 1 - V - W`, the liquid, vapour, and aqueous compositions are

.. math::
    x_i = \frac{z_i}{1 + V(K_i^V-1) + W(K_i^W-1)}, \quad y_i = K_i^V x_i, \quad w_i = K_i^W x_i

The vapour and aqueous mole fractions :math:`(V,W)` are the roots of

.. math::
    f_V = \sum_i (K_i^V-1) x_i = 0, \qquad f_W = \sum_i (K_i^W-1) x_i = 0

solved by a bounded Newton iteration. If a phase fraction goes negative or Newton fails, the solver falls back to the two-phase Rachford–Rice problem on the remaining pair (oil–gas, oil–water, or gas–water with :math:`K_i^{GW} = K_i^V / K_i^W`). Derivatives of :math:`V` and :math:`W` with respect to :math:`P`, :math:`T`, and :math:`z` follow from the implicit function theorem on the 2×2 Jacobian of :math:`(f_V, f_W)`.

Recommendations and pitfalls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The K-value flash is highly computationally efficient and is recommended for thermal recovery processes (like steam injection) or relatively stable black-oil systems where fluid compositions remain within a narrow, predictable envelope.

Pitfalls: In reality, thermodynamic K-values are heavily dependent on the overall mixture composition, not just pressure and temperature. Relying on tabulated K-values implies assuming a fixed compositional path. Care should be taken when using this model for highly variable compositional processes such as miscible gas injection (e.g., CO2 or hydrocarbon gas flooding), or strong depletion of gas condensates, where composition dependency is critical to accurate phase behaviour.

Model parameters
~~~~~~~~~~~~~~~~

The K-value flash fluid models are assigned using specific XML blocks that pair the flash model with viscosity and density models.

* ``kValueTables``: List of function names linking to the pre-tabulated K-values. For an :math:`N`-phase system, this requires :math:`(N-1) \times N_c` function names, ordered by secondary phase then component (all vapour–liquid K, then all aqueous–liquid K).
* ``pressureCoordinates``: List of pressure values for interpolation grid generation (optional). Unit: [Pa].
* ``temperatureCoordinates``: List of temperature values for interpolation grid generation (optional). Unit: [K].

Thermal K-value fluids
~~~~~~~~~~~~~~~~~~~~~~

Isothermal catalog names remain ``Compositional{Flash}Fluid{Viscosity}``. Thermal variants insert ``Thermal`` before ``Fluid`` and attach the compositional enthalpy model:

* ``CompositionalTwoPhaseKValueThermalFluidLohrenzBrayClark``
* ``CompositionalTwoPhaseKValueThermalFluidPhillipsBrine``
* ``CompositionalThreePhaseKValueThermalFluidLohrenzBrayClark``

These require ``componentHeatCapacityCoefficients`` as documented in the enthalpy model. Do not pair three-phase K-value flash with immiscible-water density if CO2 must dissolve in water; aqueous density should be compositional (or Phillips brine) so dissolved CO2 affects density.