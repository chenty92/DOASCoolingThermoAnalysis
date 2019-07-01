within ;
package HVACThermo
  model COP "COP Calculator of the HVAC system"

    // Import the package
    package Medium = Modelica.Media.Air.MoistAir;
    package Unit = Modelica.SIunits;

    // Inputs
    Modelica.Blocks.Interfaces.RealInput q_c
      "the total cooling load of the HVAC system in kW"
      annotation (Placement(transformation(extent={{-140,40},{-100,80}})));
    Modelica.Blocks.Interfaces.RealInput w_c
      "the total work done by the HVAC system in kW"
      annotation (Placement(transformation(extent={{-140,-80},{-100,-40}})));

    // Outputs
    Modelica.Blocks.Interfaces.RealOutput COP = eta
      "COP of the HVAC system"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));

    parameter Unit.Power deltaQ = 1.0
      "threshold control on total cooling load";
    parameter Unit.Power deltaWc = 1.0
      "threshold control on power input";
    parameter Unit.Power W_fan = 10.0 "fan power";
    Unit.Power Q = q_c * 1e3;
    Unit.Power W = w_c * 1e3;
    Unit.Efficiency eta "COP of the HVAC system";

  equation
    // Calculate COP of the HVAC system
    if (Q <= deltaQ) then
      eta = 0.0;
    else
      eta = Q / (W + W_fan);
    end if;

    annotation (Icon(graphics={
          Rectangle(extent={{-100,100},{100,-98}}, lineColor={28,108,200}),
          Ellipse(
            extent={{-48,46},{48,-50}},
            lineColor={0,0,0},
            fillColor={0,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-20,42},{-20,-46},{48,-6},{-20,42}},
            color={0,0,0},
            thickness=0.5),
          Text(
            extent={{-14,2},{24,-6}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,255,255},
            fillPattern=FillPattern.Solid,
            textString="COP"),
        Text(
          extent={{-150,140},{158,102}},
          lineColor={28,108,200},
          textString="%name")}));
  end COP;

  model DOASCoolingLoad
    "Calculate the required cooling load and the humidity level of the DOAS system"

    // Import the package
    package Medium = Modelica.Media.Air.MoistAir;
    package Unit = Modelica.SIunits;
    package Bldg = Buildings.Utilities.Psychrometrics.Functions;


    // Inputs
    Modelica.Blocks.Interfaces.RealInput q_lat(min=0)
      "indoor latent heat in J"
      annotation (Placement(transformation(extent={{-140,70},{-100,110}})));
    Modelica.Blocks.Interfaces.RealInput t_1(min=0)
      "outdoor air temperature in degree C"
      annotation (Placement(transformation(extent={{-140,26},{-100,66}})));
    Modelica.Blocks.Interfaces.RealInput t_3(min=0)
      "indoor air temperature in degree C"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-66},{-100,-26}})));
    Modelica.Blocks.Interfaces.RealInput rh_3(min=0)
      "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-110},{-100,-70}})));

    // Outputs
    Modelica.Blocks.Interfaces.RealOutput q_c = Q_c / 1e3
      "cooling loads for the DOAS system in kW"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput w_set=
      if w_2 >= w_2min then w_2 else w_2min
      "supply air humidity ratio"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Blocks.Interfaces.RealOutput deltaw_set=
      if w_2 >= w_2min then w_1 - w_2 else w_1 - w_2min
      "humidity removal target"
      annotation (Placement(transformation(extent={{100,-70},{120,-50}})));

    // Parameter of the model
    parameter Unit.Pressure P_1 = 101325 "outdoor air pressure";
    parameter Unit.Pressure P_3 = 101325 "indoor air pressure";
    parameter Unit.VolumeFlowRate V_a = 1.0
      "volumetric flow rate of the moist air";
    parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
    parameter Real phi_min = 0.1 "humidity control on supply air";

    // Conversion of inputs to physical parameters
    Unit.Power Q_lat = q_lat / 3.6e3
      "latent heat generated in the indoor environment";
    Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
      "outdoor air temperature";
    Real phi_1 = rh_1 / 100  "outdoor air relative humidity";
    Unit.Temperature T_3(displayUnit = "degC") = Unit.Conversions.from_degC(t_3)
      "indoor air temperature";
    Real phi_3 = rh_3 / 100 "indoor air relative humidity";

    // Flow rate of the moist air
    Unit.Density rho_a = Medium.density(State_1) "density of the moist air";
    Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";

    // State 1: Outdoor air
    Real X_1 = Medium.massFraction_pTphi(P_1, T_1, phi_1)
      "mass fraction of water vapor in the outdoor air";
    Medium.ThermodynamicState State_1 = Medium.setState_pTX(P_1, T_1, {X_1})
      "thermodynamic state of outdoor air";
    Real w_1 = X_1 / (1 - X_1) "outdoor air humidity ratio";
    Unit.SpecificEnthalpy h_1 = Medium.specificEnthalpy(State_1)
      "outdoor air specific enthalpy";

    // State 2: Supply air
    Unit.Temperature T_2(displayUnit = "degC") = T_3 "supply air temperature";
    Real w_2(min = 0) "humidity ratio of supply air";
    Real X_2 = w_2 / (1 + w_2) "mass fraction of water vapor in the supply air";
    Medium.ThermodynamicState State_2 = Medium.setState_pTX(P_3, T_2, {X_2})
      "thermodynamic state of supply air";
    Unit.SpecificEnthalpy h_2 = Medium.specificEnthalpy(State_2)
      "supply air specific enthalpy";

    // State 3: indoor air
    Real X_3 = Medium.massFraction_pTphi(P_3, T_3, phi_3)
      "mass fraction of water vapor in the indoor air";
    Medium.ThermodynamicState State_3 = Medium.setState_pTX(P_3, T_3, {X_3})
      "thermodynamic state of indoor air";
    Real w_3 = X_3 / (1 - X_3) "outdoor air humidity ratio";
    Unit.SpecificEnthalpy h_3 = Medium.specificEnthalpy(State_3)
      "indoor air specific enthalpy";
    Unit.SpecificEnthalpy h_wfg = Medium.enthalpyOfVaporization(T_3)
        "condensation heat of water vapor";

    // The total cooling load on the DOAS system
    Unit.Power Q_c "the total cooling load on the DOAS system";

    // Define the upper limit on latent load removal & total cooling load
    Unit.Power Q_cmax "maximum latent load removal rate";
    Real X_2min = Medium.massFraction_pTphi(P_3, T_3, phi_min)
      "minimum mass fraction of water vapor in the supply air";
    Real w_2min = X_2min / (1 - X_2min) "minimum humidity ratio in supply air";
    Medium.ThermodynamicState State_2min = Medium.setState_pTX(P_3, T_2, {X_2min})
      "thermodynamic state of supply air with minimum humidity ratio";
    Unit.SpecificEnthalpy h_2min = Medium.specificEnthalpy(State_2min)
      "specific enthalpy of supply air with minimum humidity ratio";


  equation
    // Remove the latent cooling load for the indoor environment
    Q_lat = m_a * (w_3 - w_2) * h_wfg;

    // The overall cooling load on the HVAC system
    Q_cmax = max(deltaQ, m_a * (h_1 - h_2min));
    Q_c = min(Q_cmax, max(deltaQ, m_a * (h_1 - h_2)));

    annotation (Documentation(info="<html>
<h4>Cooling Load Calculator for the Dedicated Outdoor Air System (DOAS)</h4>
<p><span style=\"font-family: Arial;\">The calculator consumes in the outdoor and target indoor air psychrometric conditions, and exhibits the target humidity ratio removal rate and the overall cooling load for the DOAS system.</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"),   Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Text(
            extent={{-12,78},{8,54}},
            lineColor={238,46,47},
            fillColor={0,0,0},
            fillPattern=FillPattern.None,
            textString="Q"),
          Rectangle(
            extent={{-100,2},{100,-2}},
            lineColor={102,44,145},
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-106,6},{-94,-6}},
            lineColor={238,46,47},
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{94,6},{106,-6}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-10,-52},{10,-76}},
            lineColor={238,46,47},
            fillColor={0,0,0},
            fillPattern=FillPattern.None,
            textString="w"),
        Text(
          extent={{-124,140},{134,104}},
          lineColor={28,108,200},
          textString="%name")}));
  end DOASCoolingLoad;

  model MinWork4DOASCooling

      // Import the package
      package Medium1 = Modelica.Media.Air.MoistAir;
      package Medium2 = Modelica.Media.Water.IF97_Utilities; // Steam
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput w_set(min=0)
        "supply air humidity ratio"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "process air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealInput t_4(min=0)
        "exhaust air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
      Modelica.Blocks.Interfaces.RealInput rh_4(min=0)
        "exhaust air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput w_min_l = W_min_l / 1e3
        "the minimum work of the DOAS cooling system in kW if liquid water"
        annotation (Placement(transformation(extent={{100,50},{120,70}}),
          iconTransformation(extent={{100,50},{120,70}})));
      Modelica.Blocks.Interfaces.RealOutput w_min_sa = W_min_sa / 1e3
        "the minimum work of the DOAS cooling system in kW if saturated air"
        annotation (Placement(transformation(extent={{100,10},{120,30}}),
          iconTransformation(extent={{100,10},{120,30}})));
      Modelica.Blocks.Interfaces.RealOutput w_min_erw = W_min_erw / 1e3
      "the minimum work of the DOAS cooling system in kW if ERW is implemented"
        annotation (Placement(transformation(extent={{100,-70},{120,-50}}),
          iconTransformation(extent={{100,-70},{120,-50}})));
      Modelica.Blocks.Interfaces.RealOutput w_min_hw = W_min_hw / 1e3
      "the minimum work of the DOAS cooling system in kW if HW is implemented"
        annotation (Placement(transformation(extent={{100,-30},{120,-10}}),
          iconTransformation(extent={{100,-30},{120,-10}})));


      // Parameter of the model
      parameter Unit.Pressure P_1 = 101325 "outdoor air pressure";
      parameter Unit.Pressure P_2 = 101325 "indoor air pressure";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.Efficiency epsilon_s = 0.80 "sensible heat recovery ratio";
      parameter Unit.Efficiency epsilon_l = 0.77 "latent heat recovery ratio";
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "outdoor air temperature";
      Real phi_1 = rh_1 / 100  "outdoor air relative humidity";
      Unit.Temperature T_4(displayUnit = "degC") = Unit.Conversions.from_degC(t_4)
        "indoor air temperature";
      Real phi_4 = rh_4 / 100 "indoor air relative humidity";

      // Flow rate of the moist air
      Unit.Density rho_a = Medium1.density(State_1) "density of the moist air";
      Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";

      // State 1: Outdoor air
      Real X_1 = Medium1.massFraction_pTphi(P_1, T_1, phi_1)
        "mass fraction of water vapor in the outdoor air";
      Medium1.ThermodynamicState State_1 = Medium1.setState_pTX(P_1, T_1, {X_1})
         "thermodynamic state of outdoor air";
      Real w_1 = X_1 / (1 - X_1) "outdoor air humidity ratio";
      Unit.SpecificEnergy ksi_1 "outdoor air specific exergy";
      Unit.Pressure P_1v = Bldg.pW_X(X_1) "water vapor pressure at State 1";

      // State 2: Supply air
      Unit.Temperature T_2(displayUnit = "degC") = T_4 "supply air temperature";
      Real w_2(min = 0) = w_set "humidity ratio of supply air";
      Real X_2 = w_2 / (1 + w_2) "mass fraction of water vapor in the supply air";
      Unit.SpecificEnergy ksi_2 "supply air specific exergy";

      // State 3: Liquid water
      Unit.Pressure P_3 = P_1 "liquid water pressure";
      Unit.Temperature T_3(displayUnit = "degC") = T_1
        "liquid water temperature";
      Unit.SpecificEnthalpy h_3 = Medium2.h_pT(P_3, T_3)
        "water vapor specific enthalpy";
      Unit.SpecificEntropy s_3 = Medium2.s_pT(P_3, T_3)
        "water vapor specific entropy";
      Unit.SpecificEnergy ksi_3 "liquid water specific exergy";

      // State 3e: Environment thermodynamic state for the water
      Unit.Pressure P_3e = Bldg.pW_X(X_1) "water vapor pressure";
      Unit.SpecificEnthalpy h_3e = Medium2.h_pT(P_3e, T_3)
        "water vapor specific enthalpy";
      Unit.SpecificEntropy s_3e = Medium2.s_pT(P_3e, T_3)
        "water vapor specific entropy";

      // State 4: indoor air
      Real X_4 = Medium1.massFraction_pTphi(P_2, T_4, phi_4)
        "mass fraction of water vapor in the indoor air";
      Medium1.ThermodynamicState State_4 = Medium1.setState_pTX(P_2, T_4, {X_4})
         "thermodynamic state of outdoor air";
      Real w_4 = X_4 / (1 - X_4) "indoor air humidity ratio";
      Unit.SpecificEnergy ksi_4 "indoor air specific exergy";
      Unit.Pressure P_4v = Bldg.pW_X(X_4) "water vapor pressure at State 4";
      Real X_4max = Medium1.Xsaturation(State_4)
        "mass fraction of saturated indoor air";
      Real w_4max = X_4max / (1 - X_4max) "saturated indoor air humidity ratio";


      // State 5: exhaust air for enthalpy recovery wheel
      Unit.Temperature T_5(displayUnit = "degC") "exhaust air temperature";
      Real w_5(min = 0) "exhaust air humidity ratio";
      Unit.SpecificEnergy ksi_5 "exhaust air specific exergy";

      // State 6: exhaust air for other systems
      Unit.Temperature T_6(displayUnit = "degC") = T_4 "exhaust air temperature";
      Real w_6 = min(w_4 + w_sep, w_4max) "exhaust air humidity ratio";
      Unit.SpecificEnergy ksi_6 "exhaust air specific exergy for saturated air";

      // State 7: exhaust air for heat recovery wheel
      Unit.Temperature T_7(displayUnit = "degC") "exhaust air temperature";
      Real w_7 = w_4 "exhaust air humidity ratio";
      Unit.SpecificEnergy ksi_7 "exhaust air specific exergy";


      // Minimum work input
      Real w_sep(min = 0) = w_1 - w_2 "seperation humidity ratio";
      Unit.Power W_min_sa "minimum work input for saturated air";
      Unit.Power W_min_l "minimum work input for liquid water";
      Unit.Power W_min_erw "minimum work input for exhaust air after ERW";
      Unit.Power W_min_hw "minimum work input for exhaust air after ERW";

  equation

      // Definition of effectiveness
      if (T_1 <= T_4 + deltaT or w_1 <= w_4 + deltaw or P_1v <= P_4v + deltaP) then
        T_5 = T_4;
        w_5 = w_4;
      else
        T_5 = T_4 + (T_1 - T_4) * epsilon_s;
        w_5 = w_4 + (w_1 - w_4) * epsilon_l;
      end if;

      if (T_1 <= T_4 + deltaT) then
        T_7 = T_4;
      else
        T_7 = T_4 + (T_1 - T_4) * epsilon_s;
      end if;

      // exergy
      ksi_1 = exergy_moistAir(T_1, w_1, T_1, w_1);
      ksi_2 = exergy_moistAir(T_2, w_2, T_1, w_1);
      ksi_3 = (h_3 - h_3e) - T_1 * (s_3 - s_3e);
      ksi_4 = exergy_moistAir(T_4, w_4, T_1, w_1);
      ksi_5 = exergy_moistAir(T_5, w_5, T_1, w_1);
      ksi_6 = exergy_moistAir(T_6, w_6, T_1, w_1);
      ksi_7 = exergy_moistAir(T_7, w_7, T_1, w_1);

      // Minimum work input
      W_min_l = m_a * ksi_2 + m_a * w_sep * ksi_3 - m_a * ksi_1;
      W_min_erw = m_a * ksi_2 + m_a * ksi_5 + m_a * (w_sep - (w_5 - w_4)) * ksi_3
              - m_a * ksi_1 - m_a * ksi_4;
      W_min_sa = m_a * ksi_2 + m_a * ksi_6 - m_a * ksi_1 - m_a * ksi_4;
      W_min_hw = m_a * ksi_2 + m_a * ksi_7 + m_a * w_sep * ksi_3
              - m_a * ksi_1 - m_a * ksi_4;

    annotation (Documentation(info="<html>
<h4>Cooling Load Calculator for the Dedicated Outdoor Air System (DOAS)</h4>
<p><span style=\"font-family: Arial;\">The calculator consumes in the outdoor and target indoor air psychrometric conditions, and exhibits the target humidity ratio removal rate and the overall cooling load for the DOAS system.</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"),   Icon(graphics={
          Rectangle(
            extent={{0,-20},{100,-24}},
            lineColor={102,44,145},
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-100,4},{0,-4}},
            lineColor={102,44,145},
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{0,22},{100,18}},
            lineColor={102,44,145},
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-106,6},{-94,-6}},
            lineColor={238,46,47},
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{94,26},{106,14}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-10,-52},{10,-76}},
            lineColor={238,46,47},
            fillColor={0,0,0},
            fillPattern=FillPattern.None,
            textString="w"),
        Text(
          extent={{-124,140},{134,104}},
          lineColor={28,108,200},
          textString="%name"),
          Ellipse(
            extent={{94,-16},{106,-28}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid)}));
  end MinWork4DOASCooling;

  model TWetBulb

    // Import the package
    package Medium = Modelica.Media.Air.MoistAir;
    package Unit = Modelica.SIunits;

    // Inputs
    Modelica.Blocks.Interfaces.RealInput t(min=0)
      "outdoor air temperature in degree C"
      annotation (Placement(transformation(extent={{-140,30},{-100,70}})));
    Modelica.Blocks.Interfaces.RealInput rh(min=0)
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-70},{-100,-30}})));

    // Output
    Modelica.Blocks.Interfaces.RealOutput t_wb = Unit.Conversions.to_degC(T_wb)
      "wetbulb temperature in degree C"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));

    // Parameter
    Unit.Pressure P = 101325 "Ambient pressure in Pa";

    // Conversions
    Unit.Temperature T(displayUnit = "degC") = Unit.Conversions.from_degC(t)
      "outdoor air temperature";
    Real phi = rh / 100 "outdoor air relative humidity";

    // Calculate wet bulb temperature
    Buildings.Utilities.Psychrometrics.TWetBul_TDryBulPhi wetBul(Medium);
    Unit.Temperature TDryBul "dry-bulb temperature";
    Unit.Pressure p "pressure";
    Real relHum "relative humidity";
    Unit.Temperature T_wb(displayUnit = "degC")= wetBul.TWetBul;

  equation
    // Calculate wet bulb temperature
    TDryBul = T;
    relHum = phi;
    p = P;
    TDryBul = wetBul.TDryBul;
    p = wetBul.p;
    relHum = wetBul.phi;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
        Text(
          extent={{-140,140},{148,104}},
          lineColor={28,108,200},
          textString="%name"),
          Rectangle(extent={{-100,100},{100,-100}},lineColor={28,108,200}),
          Ellipse(
            extent={{-60,60},{60,-60}},
            lineColor={28,108,200},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{10,-4},{-12,4}},
            lineColor={255,255,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.None,
            textString="T + RH -> Twb")}),                         Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end TWetBulb;

  model EnthalpyRecoveryWheel "ERW Component"

      // Import the package
      package Medium = Modelica.Media.Air.MoistAir;
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "hot air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_2(min=0)
        "cool air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,8},{-100,48}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "hot air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-48},{-100,-8}})));
      Modelica.Blocks.Interfaces.RealInput rh_2(min=0)
        "cool air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput t_3 = T_3 - 273.15
        "dehumidified and cooled air temperature"
        annotation (Placement(transformation(extent={{100,36},{120,56}}),
          iconTransformation(extent={{100,36},{120,56}})));
      Modelica.Blocks.Interfaces.RealOutput t_4 = T_4 - 273.15
        "humidified exhaust air temperature in degree C"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
          iconTransformation(extent={{100,-10},{120,10}})));
      Modelica.Blocks.Interfaces.RealOutput rh_3 = phi_3 * 100
        "dehumidified and cooled air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-54},{120,-34}}),
          iconTransformation(extent={{100,-54},{120,-34}})));
      Modelica.Blocks.Interfaces.RealOutput rh_4 = phi_4 * 100
        "humidified exhaust air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-100},{120,-80}}),
          iconTransformation(extent={{100,-100},{120,-80}})));
      Modelica.Blocks.Interfaces.RealOutput w_p = W_p / 1e3
        "fan power of the ERV system in kW"
        annotation (Placement(transformation(extent={{100,80},{120,100}}),
          iconTransformation(extent={{100,80},{120,100}})));


      // Parameter of the model
      parameter Unit.Pressure P = 101325 "air pressure";
      parameter Unit.Efficiency epsilon_s(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Efficiency epsilon_l(max = 1.0) = 0.75
        "effectiveness of latent heat recovery";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in ERW to supply air";
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";
      parameter Unit.Pressure P_fan = 100
        "pressure loss across the ERV system";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "hot air temperature";
      Real phi_1 = rh_1 / 100 "hot air relative humidity";
      Unit.Temperature T_2(displayUnit = "degC") = Unit.Conversions.from_degC(t_2)
        "cool air temperature";
      Real phi_2 = rh_2 / 100 "cool air relative humidity";

      // State 1: Hot air stream
      Real X_1 = Medium.massFraction_pTphi(P, T_1, phi_1)
        "mass fraction of water vapor in the hot air";
      Medium.ThermodynamicState State_1 = Medium.setState_pTX(P, T_1, {X_1})
        "thermodynamic state of hot air";
      Real w_1 = X_1 / (1 - X_1) "hot air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium.specificEnthalpy(State_1)
        "hot air specific enthalpy";
      Unit.Pressure P_1sat = Medium.saturationPressure(T_1)
        "outdoor air satuation pressure";
      Unit.Pressure P_1v "water vapor pressure at State 1";

      // State 2: Cool air stream
      Real X_2 = Medium.massFraction_pTphi(P, T_2, phi_2)
        "mass fraction of water vapor in the cool air";
      Medium.ThermodynamicState State_2 = Medium.setState_pTX(P, T_2, {X_2})
        "thermodynamic state of cool air";
      Real w_2 = X_2 / (1 - X_2) "cool air humidity ratio";
      Unit.SpecificEnthalpy h_2 = Medium.specificEnthalpy(State_2)
        "cool air specific enthalpy";
      Unit.Pressure P_2sat = Medium.saturationPressure(T_2)
        "outdoor air satuation pressure";
      Unit.Pressure P_2v "water vapor pressure at State 2";

      // State 3: Dehumidified and cooled air stream
      Unit.Temperature T_3(displayUnit = "degC")
        "dehumidified and cooled air temperature";
      Real w_3 "dehumidified and cooled air humidity ratio";
      Real X_3 = w_3 / (1 + w_3) "mass fraction of water vapor in the cool air";
      Medium.ThermodynamicState State_3 = Medium.setState_pTX(P, T_3, {X_3})
        "thermodynamic state of cool air";
      Unit.SpecificEnthalpy h_3 = Medium.specificEnthalpy(State_3)
        "dehumidified and cooled air specific enthalpy";
      Real phi_3 = Medium.relativeHumidity(State_3)
        "dehumidified and cooled air relative humidity";

      // State 4: hot exhaust air stream
      Real w_4 "hot exhaust air humidity ratio";
      Real X_4 = w_4 / (1 + w_4)
        "mass fraction of water vapor in the hot exhaust air";
      Unit.SpecificEnthalpy h_4 "hot exhaust air specific enthalpy";
      Medium.ThermodynamicState State_4 = Medium.setState_phX(P, h_4, {X_4})
        "thermodynamic state of hot exhaust air";
      Unit.Temperature T_4 = Medium.temperature(State_4)
        "hot exhaust air temperature";
      Real phi_4 = Medium.relativeHumidity(State_4)
        "hot exhaust air relative humidity";

      // Fan Power output
      Unit.Power W_p "fan power in W";

  equation
      // Water vapor pressure at State 1
      phi_1 = P_1v / P_1sat;
      phi_2 = P_2v / P_2sat;

      // Definition of effectiveness
      // epsilon = m_a * Cp_a * (Torw_1 - Torw_3) /
      //          (m_a * r * Cp_a * (Torw_1 - Torw_2))
      if (T_1 <= T_2 + deltaT or w_1 <= w_2 + deltaw or P_1v <= P_2v + deltaP) then
        T_3 = T_1;
        w_3 = w_1;
        W_p = 0.0;
      else
        epsilon_s = (T_1 - T_3) / r / (T_1 - T_2);
        epsilon_l = (w_1 - w_3) / r / (w_1 - w_2);
        W_p = P_fan * V_a;
      end if;

      // Mass conservation for the water content
      w_1 + r * w_2 = w_3 + r * w_4;

      // Energy conservation for the whole system
      h_1 + r * h_2 = h_3 + r * h_4;

    annotation (Line(points={{110,30},{108,30},{108,30},{110,30}},
          color={0,0,127}),
      Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">ERW (Enthalpy Recovery Wheel) Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the ERW given the effectiveness of the sensible and latent heat recovery. The model also includes the threshold control on the ERW, i.e. when the temperature of inlet hot air is smaller than the threshold above the inlet cool air, the sensible heat recovery system is off; whnn the humidity ratio differerence or water vapor patial pressure of the inlet humid air is smaller than the threshold above the inlet dry air, the latent heat recovery system is off.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Cool air inlet; State 3: Hot air outlet (dehumidifed and cooled air); State 4: Cool air outlet (humidified air)</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"),
      Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-40,80},{40,-80}},
            lineColor={28,108,200},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={238,46,47}),
          Ellipse(
            extent={{-100,60},{-82,42}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={238,46,47}),
          Ellipse(
            extent={{82,-40},{100,-58}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={127,12,91}),
          Ellipse(
            extent={{82,60},{100,42}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={28,108,200}),
          Ellipse(
            extent={{-100,-40},{-82,-58}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={102,44,145}),
          Line(
            points={{-82,50},{-40,50},{40,-50},{82,-50}},
            color={28,108,200},
            thickness=1),
          Line(
            points={{82,50},{40,50},{-40,-50},{-82,-50}},
            color={162,29,33},
            thickness=1,
            origin={0,0},
            rotation=360),
        Text(
          extent={{-122,138},{132,102}},
          lineColor={28,108,200},
          textString="%name")}));
  end EnthalpyRecoveryWheel;

  model HeatRecoveryWheel "Heat Recovery Wheel"

      // Import the package
      package Medium = Modelica.Media.Air.MoistAir;
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "hot air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_2(min=0)
        "cool air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,8},{-100,48}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "hot air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-48},{-100,-8}})));
      Modelica.Blocks.Interfaces.RealInput rh_2(min=0)
        "cool air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput t_3 = T_3 - 273.15
        "cooled process air temperature"
        annotation (Placement(transformation(extent={{100,34},{120,54}}),
          iconTransformation(extent={{100,34},{120,54}})));
      Modelica.Blocks.Interfaces.RealOutput t_4 = T_4 - 273.15
        "hot exhaust air temperature in degree C"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
          iconTransformation(extent={{100,-10},{120,10}})));
      Modelica.Blocks.Interfaces.RealOutput rh_3 = phi_3 * 100
        "cooled process air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-54},{120,-34}}),
          iconTransformation(extent={{100,-54},{120,-34}})));
      Modelica.Blocks.Interfaces.RealOutput rh_4 = phi_4 * 100
        "hot exhaust air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-100},{120,-80}}),
          iconTransformation(extent={{100,-100},{120,-80}})));
      Modelica.Blocks.Interfaces.RealOutput w_p = W_p / 1e3
        "fan power of the HRW system in kW"
        annotation (Placement(transformation(extent={{100,78},{120,98}}),
          iconTransformation(extent={{100,78},{120,98}})));



      // Parameter of the model
      parameter Unit.Pressure P = 101325 "air pressure";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in HW to supply air";
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the HRW system";
      parameter Unit.Pressure P_fan = 100
        "pressure loss across the HRW system";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "hot air temperature";
      Real phi_1 = rh_1 / 100 "hot air relative humidity";
      Unit.Temperature T_2(displayUnit = "degC") = Unit.Conversions.from_degC(t_2)
        "cool air temperature";
      Real phi_2 = rh_2 / 100 "cool air relative humidity";

      // State 1: Hot air stream
      Real X_1 = Medium.massFraction_pTphi(P, T_1, phi_1)
        "mass fraction of water vapor in the hot air";
      Medium.ThermodynamicState State_1 = Medium.setState_pTX(P, T_1, {X_1})
        "thermodynamic state of hot air";
      Real w_1 = X_1 / (1 - X_1) "hot air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium.specificEnthalpy(State_1)
        "hot air specific enthalpy";

      // State 2: Cool air stream
      Real X_2 = Medium.massFraction_pTphi(P, T_2, phi_2)
        "mass fraction of water vapor in the cool air";
      Medium.ThermodynamicState State_2 = Medium.setState_pTX(P, T_2, {X_2})
        "thermodynamic state of cool air";
      Real w_2 = X_2 / (1 - X_2) "cool air humidity ratio";
      Unit.SpecificEnthalpy h_2 = Medium.specificEnthalpy(State_2)
        "cool air specific enthalpy";

      // State 3: Process air stream
      Unit.Temperature T_3(displayUnit = "degC")
        "process air temperature";
      Real X_3 = X_1 "mass fraction of water vapor in the cool air";
      Medium.ThermodynamicState State_3 = Medium.setState_pTX(P, T_3, {X_3})
        "thermodynamic state of cool air";
      Unit.SpecificEnthalpy h_3 = Medium.specificEnthalpy(State_3)
        "process air specific enthalpy";
      Real phi_3 = Medium.relativeHumidity(State_3)
        "process air relative humidity";

      // State 4: hot exhaust air stream
      Real X_4 = X_2 "mass fraction of water vapor in the hot exhaust air";
      Unit.SpecificEnthalpy h_4(start = 6e4) "hot exhaust air specific enthalpy";
      Medium.ThermodynamicState State_4 = Medium.setState_phX(P, h_4, {X_4})
        "thermodynamic state of hot exhaust air";
      Unit.Temperature T_4 = Medium.temperature(State_4)
        "hot exhaust air temperature";
      Real phi_4 = Medium.relativeHumidity(State_4)
        "hot exhaust air relative humidity";

      // Fan Power output
      Unit.Power W_p "fan power in W";


  equation

      // Definition of effectiveness
      // epsilon = m_a * Cp_a * (T_1 - T_3) /
      //          (m_a * r * Cp_a * (T_1 - T_2))
      if (T_1 <= T_2 + deltaT) then
        T_3 = T_1;
        W_p = 0.0;
      else
        epsilon = (T_1 - T_3) / r / (T_1 - T_2);
        W_p = P_fan * V_a;
      end if;

      // Energy conservation for the whole system
      h_1 + r * h_2 = h_3 + r * h_4;

    annotation (Line(points={{110,30},{108,30},{108,30},{110,30}},
          color={0,0,127}),
      Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">ERW (Enthalpy Recovery Wheel) Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the ERW given the effectiveness of the sensible and latent heat recovery. The model also includes the threshold control on the ERW, i.e. when the temperature of inlet hot air is smaller than the threshold above the inlet cool air, the sensible heat recovery system is off; whnn the humidity ratio differerence or water vapor patial pressure of the inlet humid air is smaller than the threshold above the inlet dry air, the latent heat recovery system is off.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Cool air inlet; State 3: Hot air outlet (dehumidifed and cooled air); State 4: Cool air outlet (humidified air)</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"),
      Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-40,80},{40,-80}},
            lineColor={28,108,200},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={238,46,47}),
          Ellipse(
            extent={{-100,60},{-82,42}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={238,46,47}),
          Ellipse(
            extent={{82,-40},{100,-58}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={127,12,91}),
          Ellipse(
            extent={{82,60},{100,42}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={28,108,200}),
          Ellipse(
            extent={{-100,-40},{-82,-58}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={102,44,145}),
          Line(
            points={{-82,50},{-40,50},{40,-50},{82,-50}},
            color={28,108,200},
            thickness=1),
          Line(
            points={{82,50},{40,50},{-40,-50},{-82,-50}},
            color={162,29,33},
            thickness=1,
            origin={0,0},
            rotation=360),
        Text(
          extent={{-122,138},{132,102}},
          lineColor={28,108,200},
          textString="%name")}));
  end HeatRecoveryWheel;

  model GCHX "Ground-Coupled Heat Exchanger"


      // Import the package
      package Medium = Modelica.Media.Air.MoistAir;
      package Medium1 = Modelica.Media.Water.WaterIF97_R1pT; // Liquid water
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;
      package Data = Buildings.Fluid.Chillers.Data.ElectricEIR;
      package Func = Buildings.Utilities.Math.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput q_c(min=0)
        "additional cooling load in kW"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
      "hot air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput t_2(min=0)
        "cool air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "hot air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
      Modelica.Blocks.Interfaces.RealInput t_wb(min=0)
      "wet bulb temperature of outdoor air in degree C"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));


      // Outputs
      Modelica.Blocks.Interfaces.RealOutput w_c = W_c / 1e3
        "the total work of the GCHX system in kW"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
          iconTransformation(extent={{100,-10},{120,10}})));

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "air pressure";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Diameter d_p = 2.5e-2 "diameter of pipe";
      parameter Unit.Height H = 30 "depth of the borehole";
      parameter Unit.Efficiency epsilon_b(max = 1.0) = 0.80
        "effectiveness of ground heat exchanger";
      parameter Real n = 1.0 "number of parallel boreholes";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the GCHX system";
      parameter Unit.Temperature t_g(min = 0) = 10.0
        "annual average ground temperature in degC";
      parameter Real gCHXisOn = 1.0
        "whether groud heat exchanger is used in the analysis, 1 represent on";
      parameter Real r_R = 0.1
        "ratio of thermal resistance of borehole to the ground";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "hot air temperature in K";
      Real phi_1 = rh_1 / 100 "hot air relative humidity";
      Unit.Temperature T_2(displayUnit = "degC") = Unit.Conversions.from_degC(t_2)
        "cool air temperature in K";
      Unit.Temperature T_g = Unit.Conversions.from_degC(t_g)
        "annual average ground temperature in K";
      Unit.Power Q_c = q_c * 1e3 "additional cooling load in W";
      Unit.Temperature T_wb(displayUnit = "degC") = Unit.Conversions.from_degC(t_wb)
        "outdoor air wet bulb temperature in K";

      // Flow rate of the moist air
      Unit.Density rho_a = Medium.density(State_1) "density of the moist air";
      Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";
      Unit.SpecificHeatCapacity Cp_a = Medium.specificHeatCapacityCp(State_1)
        "process air specific heat capacity";

      // Flow rate of the coolant
      Unit.Density rho_w = Medium1.density(State_3)
        "density of the chilled water";
      Unit.Velocity v_w(start = 1e-2, min = 0.0) "flow rate of the coolant";
      Unit.MassFlowRate m_w(start = 1e-2, min = 0.0)
        "mass flow rate of the coolant";
      Unit.SpecificHeatCapacity Cp_w = Medium1.specificHeatCapacityCp(State_3)
        "coolant specific heat capacity";
      Unit.DynamicViscosity mu_w = Medium1.dynamicViscosity(State_3)
        "viscosity of coolant";

      // State 1: Hot air stream
      Real X_1 = Medium.massFraction_pTphi(P, T_1, phi_1)
        "mass fraction of water vapor in the hot air";
      Medium.ThermodynamicState State_1 = Medium.setState_pTX(P, T_1, {X_1})
        "thermodynamic state of hot air";
      Real w_1 = X_1 / (1 - X_1) "hot air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium.specificEnthalpy(State_1)
        "hot air specific enthalpy";

      // State 2: Cool air stream
      Real w_2 = w_1 "humidity ratio in the supply air";
      Real X_2 = w_2 / (1 + w_2)
        "mass fraction of water vapor in the supply air";
      Medium.ThermodynamicState State_2 = Medium.setState_pTX(P, T_2, {X_2})
        "thermodynamic state of cool air";
      Unit.SpecificEnthalpy h_2 = Medium.specificEnthalpy(State_2)
        "cool air specific enthalpy";

      // State 3: Coolant inlet
      Unit.Temperature T_3(displayUnit = "degC", min = 273.0)
        "coolant inlet temperature";
      Real X_3 = 1.0 "mass fraction of liquid water in the chilled water";
      Medium1.ThermodynamicState State_3 = Medium1.setState_pTX(P, T_3, {X_3})
        "thermodynamic state of coolant";
      Unit.SpecificEnthalpy h_3 = Medium1.specificEnthalpy(State_3)
        "coolant inlet specific enthalpy";

      // State 4: Coolant outlet
      Real X_4 = 1.0 "mass fraction of liquid water in the chilled water";
      Unit.Temperature T_4(displayUnit = "degC", min = 273.0)
        "coolant outlet temperature";
      Medium1.ThermodynamicState State_4 = Medium1.setState_pTX(P, T_4, {X_4})
        "thermodynamic state of coolant outlet";
      Unit.SpecificEnthalpy h_4 = Medium1.specificEnthalpy(State_4)
        "coolant outlet specific enthalpy";

      // Effectiveness definition
      Real mCpmin(unit = "W/K")
        "minimum heat capacity between moist air and chilled water";
      Unit.Temperature T_b(displayUnit = "degC") "borehole wall temperature";

      // Pump Power
      Unit.Power W_c "work done by the pump";
      Unit.Pressure deltaP "pressure drop along the pipe";
      Real f "friction factor";
      Real pi = Modelica.Constants.pi;
      Real Re "Reynolds number of pipe flow";

      // Chiller model as an alternative
      // Chiller model type
      Data.ElectricEIRChiller_Trane_CVHF_2567kW_11_77COP_VSD per;
      Unit.Temperature T_lmin(displayUnit = "degC")=
        Unit.Conversions.to_degC(per.TEvaLvgMin)
        "lower limit on leaving chilled water temperature";
      Unit.Temperature T_lmax(displayUnit = "degC")=
        Unit.Conversions.to_degC(per.TEvaLvgMax)
        "upper limit on leaving chilled water temperature";
      Unit.Temperature T_emin(displayUnit = "degC")=
        Unit.Conversions.to_degC(per.TConEntMin)
        "lower limit on entering condenser water temperature";
      Unit.Temperature T_emax(displayUnit = "degC")=
        Unit.Conversions.to_degC(per.TConEntMax)
        "upper limit on entering condenser water temperature";
      Unit.Temperature Tl = min(T_lmax, max(T_lmin, T_3))
        "actual leaving chilled water temperature";
      Unit.Temperature Te = min(T_emax, max(T_emin, T_wb + deltaT_wc))
        "actual entering condenser water temperature";
      // Chiller model paramaters
      Unit.Efficiency capFunT = Func.smoothMax(
        x1 = 1E-6,
        x2 = Func.biquadratic(a = per.capFunT, x1 = Tl, x2 = Te),
        deltaX = 1E-7);
      Unit.Efficiency EIRFunT = Func.biquadratic(a = per.EIRFunT, x1 = Tl, x2 = Te);
      Unit.Efficiency EIRFunPLR = per.EIRFunPLR[1] + per.EIRFunPLR[2] * PLR +
        per.EIRFunPLR[3] * PLR^2;
      Unit.Efficiency PLR "part load ratio";
      Unit.Efficiency COP_c "substitute of cooling coil due to higher ground temperature";
      Real unit "number of parallel terminals";
      Unit.Temperature deltaT_wc(displayUnit = "degC") = 2
        "approximate constant approach temperature";


  equation

      // Definition of effectiveness
      // epsilon = m_a * Cp_a * (T_1 - T_2) /
      //          (mCpmin * (T_1 - T_3))
      // mCpmin = min(m_a * Cp_a, m_w * Cp_w) = m_a * Cp_a;
      mCpmin = m_a * Cp_a;

      // Effectiveness of heat exchanger
      epsilon = m_a * Cp_a * (T_1 - T_2) / (mCpmin * (T_1 - T_3));

      // Energy conservation for the whole system
      m_a * Cp_a * (T_1 - T_2) + Q_c = m_w * (h_4 - h_3);

      // Effectiveness of ground heat exchanger
      epsilon_b = (T_4 - T_3) / (T_4 - T_b);

      // Calculate the borehole wall temperature
      T_b = (T_3 + r_R * T_g) / (1 + r_R);

      // Pressure drop in each borehole
      deltaP = f * 1/2 * rho_w * v_w^2 * (2 * H) / d_p;
      m_w = rho_w * v_w * (n * pi * d_p^2 / 4);
      Re = rho_w * v_w * d_p / mu_w;
      if (Re <= 0) then
        f = 0;
      else
        if (Re < 2.3e3) then
          f = 64 / Re;
        else
          f = (1.82 * Modelica.Math.log10(Re) - 1.64) ^ (-2);
        end if;
      end if;

      // Pump work
      if (T_1 <= T_2 + deltaT) then
        // Inlet air must be cooler than the exhaust air
        W_c = 0.0;
        PLR = 0.0;
        COP_c = 0.0;
        unit = 0.0;
      elseif (gCHXisOn <= 1e-2 or T_g >= T_3 - deltaT) then
        // Ground temperature must be cooler than the inlet coolant water
        W_c = m_w * (h_4 - h_3) / COP_c;
        PLR = 0.5;
        PLR = m_w * (h_4 - h_3)/(-per.QEva_flow_nominal/unit)/capFunT;
        COP_c = per.COP_nominal / EIRFunT / EIRFunPLR * PLR;
      else
        W_c = deltaP * m_w / rho_w;
        PLR = 0.0;
        COP_c = 0.0;
        unit = 0.0;
      end if;

    annotation (Line(points={{110,30},{108,30},{108,30},{110,30}},
          color={0,0,127}),
      Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">ERW (Enthalpy Recovery Wheel) Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the ERW given the effectiveness of the sensible and latent heat recovery. The model also includes the threshold control on the ERW, i.e. when the temperature of inlet hot air is smaller than the threshold above the inlet cool air, the sensible heat recovery system is off; whnn the humidity ratio differerence or water vapor patial pressure of the inlet humid air is smaller than the threshold above the inlet dry air, the latent heat recovery system is off.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Cool air inlet; State 3: Hot air outlet (dehumidifed and cooled air); State 4: Cool air outlet (humidified air)</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"),
      Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-40,80},{40,-80}},
            lineColor={28,108,200},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={238,46,47}),
          Ellipse(
            extent={{-100,60},{-82,42}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={238,46,47}),
          Ellipse(
            extent={{82,-40},{100,-58}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={127,12,91}),
          Ellipse(
            extent={{82,60},{100,42}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={28,108,200}),
          Ellipse(
            extent={{-100,-40},{-82,-58}},
            lineColor={28,108,200},
            fillPattern=FillPattern.Solid,
            fillColor={102,44,145}),
          Line(
            points={{-82,50},{-40,50},{40,-50},{82,-50}},
            color={28,108,200},
            thickness=1),
          Line(
            points={{82,50},{40,50},{-40,-50},{-82,-50}},
            color={162,29,33},
            thickness=1,
            origin={0,0},
            rotation=360),
        Text(
          extent={{-122,138},{132,102}},
          lineColor={28,108,200},
          textString="%name")}));
  end GCHX;

  model Chiller "Chiller component"

    // Import the package
    package Medium = Modelica.Media.Air.MoistAir;
    package Medium1 = Modelica.Media.Water.WaterIF97_R1pT; // Liquid water
    package Unit = Modelica.SIunits;
    package Bldg = Buildings.Utilities.Psychrometrics.Functions;

    // inputs
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_4(min=0)
        "indoor air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,6},{-100,46}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "process air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-46},{-100,-6}})));
      Modelica.Blocks.Interfaces.RealInput w_set(min=0)
        "supply air humidity ratio"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

    //Outputs
    Modelica.Blocks.Interfaces.RealOutput q_c = Q_c / 1e3
      "cooling loads for chiller in kW"
      annotation (Placement(transformation(extent={{100,70},{120,90}})));
    Modelica.Blocks.Interfaces.RealOutput t_5 = Unit.Conversions.to_degC(T_5)
      "leaving chilled water temperature in degree C"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Blocks.Interfaces.RealOutput t_wc = Unit.Conversions.to_degC(T_wc)
      "entering condenser water temperature in degree C"
      annotation (Placement(transformation(extent={{100,-90},{120,-70}})));

    // Parameter of the model
    parameter Unit.Pressure P_1 = 101325 "outdoor air pressure";
    parameter Unit.Pressure P_4 = 101325 "indoor air pressure";
    parameter Unit.VolumeFlowRate V_a = 1.0
      "volumetric flow rate of the moist air";
    parameter Unit.VolumeFlowRate V_w = 1e-2
      "volumetric flow rate of the chilled water";
    parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
      "effectiveness of the heat exchanger";
    parameter Unit.Power deltaQ = 1.0 "minimum cooling load requirement";
    parameter Real t_lcmin = 1.0 "minimum leaving chilled water temperature in degC";

    // Conversion of inputs to physical parameters
    Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
      "process air temperature";
    Real phi_1 = rh_1 / 100  "process air relative humidity";
    Unit.Temperature T_4(displayUnit = "degC") = Unit.Conversions.from_degC(t_4)
      "indoor air temperature";

    // Flow rate of the moist air
    Unit.Density rho_a = Medium.density(State_1) "density of the moist air";
    Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";
    Unit.SpecificHeatCapacity Cp_a = Medium.specificHeatCapacityCp(State_1)
      "process air specific heat capacity";

    // State 1: Process air
    Real X_1 = Medium.massFraction_pTphi(P_1, T_1, phi_1)
      "mass fraction of water vapor in the process air";
    Medium.ThermodynamicState State_1 = Medium.setState_pTX(P_1, T_1, {X_1})
      "thermodynamic state of process air";
    Real w_1 = X_1 / (1 - X_1) "process air humidity ratio";
    Unit.SpecificEnthalpy h_1 = Medium.specificEnthalpy(State_1)
      "process air specific enthalpy";

    // State 2: Condensation of process air
    Real w_2 "condensation air humidity ratio";
    Real X_2 "mass fraction of water vapor in the condensation air";
    Unit.Pressure P_2 = 101325 "condensation air pressure";
    // Calculate the dewpoint temperature of process air
    Unit.Pressure pW_1 = Bldg.pW_X(X_1) "water vapor pressure of process air";
    Unit.Temperature T_2(displayUnit = "degC") = Bldg.TDewPoi_pW(pW_1)
      "Dewpoint temperature of process air";
    Medium.ThermodynamicState State_2 = Medium.setState_pTX(P_2, T_2, {X_2})
      "thermodynamic state of process air";
    Unit.SpecificEnthalpy h_2 = Medium.specificEnthalpy(State_2)
      "condensation air specific enthalpy";

    // State 3: Dehumidified cool air
    Real w_3 "dehumidified air humidity ratio";
    Real X_3 = w_3 / (1 + w_3)
      "mass fraction of water vapor in the dehumidified air";
    Unit.SpecificEnthalpy h_wfg = Medium.enthalpyOfVaporization(T_2)
      "condensation heat of water vapor";
    Unit.Pressure P_3 = 101325 "dehumidified air pressure";
    Unit.Pressure pW_3 = Bldg.pW_X(X_3)
      "water vapor pressure of dehumidified air";
    Unit.Temperature T_3(displayUnit = "degC") = Bldg.TDewPoi_pW(pW_3)
      "dewpoint temperature of dehumidified air";
    Medium.ThermodynamicState State_3 = Medium.setState_pTX(P_3, T_3, {X_3})
      "thermodynamic state of dehumidified air";
    Unit.SpecificEnthalpy h_3 = Medium.specificEnthalpy(State_3)
      "dehumidified air specific enthalpy";

    // The total cooling load on the chiller
    Unit.Power Q_c "the total cooling load on the chiller";

    // Leaving chilled water temperature
    Unit.Temperature T_5(displayUnit = "degC", start = 283.0, min = 273.0,
      max = 373.0) "chilled water temperature";
    Real X_5 = 1.0 "mass fraction of liquid water in the chilled water";
    Unit.Pressure P_5 = 101325 "chilled water pressure";
    Medium1.ThermodynamicState State_5 = Medium1.setState_pTX(P_5, T_5, {X_5});
    Unit.SpecificEnthalpy h_5 = Medium1.specificEnthalpy(State_5)
      "chilled water specific enthalpy";
    Unit.Density rho_w = Medium1.density(State_5)
      "density of the chilled water";
    Unit.MassFlowRate m_w = rho_w * V_w "mass flow rate of the moist air";
    Unit.SpecificHeatCapacity Cp_w = Medium1.specificHeatCapacityCp(State_5)
      "chilled water specific heat capacity";

    // Effectiveness definition
    Real mCpmin(unit = "W/K")
      "minimum heat capacity between moist air and chilled water";

    // Entering chilled water temperature (heated water at the exit of HX)
    Unit.Pressure P_6 = 101325 "chilled water pressure";
    Unit.SpecificEnthalpy h_6 "specific enthalpy of the heated water";
    Unit.Temperature T_6 = Medium1.temperature_ph(P_6, h_6)
      "heated water temperature";

    // Entering condenser water temperature (from cooling tower)
    Buildings.Utilities.Psychrometrics.TWetBul_TDryBulPhi wetBul(Medium);
    Unit.Temperature TDryBul "dry-bulb temperature";
    Unit.Pressure p "pressure";
    Real phi "relative humidity";
    Unit.Temperature deltaT_wc(displayUnit = "degC") = 2
      "approximate constant approach temperature";
    Unit.Temperature T_wc(displayUnit = "degC")= wetBul.TWetBul + deltaT_wc;

    // Determine the maximum cooling load state
    Unit.Temperature T_5min(displayUnit = "degC")=
        Unit.Conversions.from_degC(t_lcmin)
      "minimum available chilled water temperature";
    Unit.Temperature T_3min(displayUnit = "degC")
      "minimum dew point temperature with maximum cooling load";
    Unit.Pressure pW_3min = Bldg.pW_TDewPoi(T_3min)
      "water vapor pressure in the dehumidified air with maximum cooling load";
    Real X_3min = Bldg.X_pW(pW_3min)
      "mass fraction of water vapor in the dehumidified air with maximum cooling load";
    Real w_3min = X_3min / (1 - X_3min)
      "dehumidified air humidity ratio with maximum cooling load";


  equation
    // From State 1 to State 2, cool down the process air to condensation
    w_2 = w_1;
    X_2 = X_1;

    // The overall cooling load on the chiller
    Q_c = max(deltaQ, m_a * (h_1 - h_3));

    // Calculate the leaving chilled water temperature from HX constant
    // effectiveness model
    // mCpmin = min(m_a * Cp_a, m_w * Cp_w) = m_a * Cp_a;
    mCpmin = m_a * Cp_a;
    epsilon = m_a * Cp_a * (T_1 - T_3) / (mCpmin * (T_1 - T_5));

    // Calculate the maximum cooling load on the chiller
    epsilon = m_a * Cp_a * (T_1 - T_3min) / (mCpmin * (T_1 - T_5min));
    w_3 = max(w_3min, w_set);

    // Heat balance of the HX, to determine the State of heated water
    m_a * h_1 + m_w * h_5 = m_a * h_3 + m_w * h_6;

    // Calculate Entering condenser water temperature
    TDryBul = T_1;
    phi = phi_1;
    p = P_1;
    TDryBul = wetBul.TDryBul;
    p = wetBul.p;
    phi = wetBul.phi;

    annotation (Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Chiller Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the chiller given the outdoor air and indoor air conditions from E+. The inputs are prepared by E+ exported results; and the output of the model is used to evaluate the COP of the chiller for a certain time period of the year at a specific location.</span></p>
<p><img src=\"modelica://HVACThermo/Chiller.png\"/></p>
<p><span style=\"font-family: Arial;\">State 1: Outdoor hot air; State 2: Saturated air for condensation; State 3: Dehumidified supply air; State 4: Indoor cool air; State 5: Chiller water leaving evaporator </span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"),
      Icon(graphics={
          Rectangle(
            extent={{-100,2},{100,-2}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={217,67,180},
            fillPattern=FillPattern.HorizontalCylinder),
          Rectangle(
            extent={{-60,40},{60,-40}},
            lineColor={28,108,200},
            lineThickness=0.5,
            fillColor={30,192,200},
            fillPattern=FillPattern.Sphere),
          Line(
            points={{4,40},{-6,30},{4,20},{-6,10},{4,0},{-6,-10},{4,-20},{-6,-30},
                {4,-40}},
            color={238,46,47},
            thickness=1,
            smooth=Smooth.Bezier),
          Ellipse(
            extent={{-100,6},{-88,-6}},
            lineColor={238,46,47},
            lineThickness=0.5,
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{88,6},{100,-6}},
            lineColor={28,108,200},
            lineThickness=0.5,
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Rectangle(extent={{-100,100},{100,-100}},
                                                lineColor={28,108,200}),
        Text(
          extent={{-140,140},{148,104}},
          lineColor={28,108,200},
          textString="%name")}));
  end Chiller;

  model ChillerCOPFittingCurve "COP fitting curve function for Chiller"

    // Import the package
    package Unit = Modelica.SIunits;
    package Data = Buildings.Fluid.Chillers.Data.ElectricEIR;
    package Func = Buildings.Utilities.Math.Functions;

    // inputs
    Modelica.Blocks.Interfaces.RealInput q_c(min=0)
      "cooling loads for chiller in kW"
        annotation (Placement(transformation(extent={{-140,50},{-100,90}})));
    Modelica.Blocks.Interfaces.RealInput t_l(min=0)
      "leaving chilled water temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput t_e(min=0)
      "entering condenser water temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-90},{-100,-50}})));

    // Output
    Modelica.Blocks.Interfaces.RealOutput COP = eta
      "COP of the chiller"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.RealOutput w_c = w_tot / 1e3
      "work done by the chiller in kW"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}})));

    // Parameter of the model
    parameter Unit.Power deltaQ = 1.0 "threshold control on the cooling load";
    parameter Unit.Efficiency deltaCOP = 1e-2;

    // Chiller model type
    Data.ElectricEIRChiller_Trane_CVHF_2567kW_11_77COP_VSD per;

    // Conversion of inputs to physical parameters
    Unit.Power Q = q_c * 1e3 "the cooling load in W";
    Unit.Temperature T_lmin(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TEvaLvgMin)
      "lower limit on leaving chilled water temperature";
    Unit.Temperature T_lmax(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TEvaLvgMax)
      "upper limit on leaving chilled water temperature";
    Unit.Temperature T_emin(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TConEntMin)
      "lower limit on entering condenser water temperature";
    Unit.Temperature T_emax(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TConEntMax)
      "upper limit on entering condenser water temperature";

    // Regression Model for the Chiller
    Unit.Temperature Tl = min(T_lmax, max(T_lmin, t_l))
      "actual leaving chilled water temperature";
    Unit.Temperature Te = min(T_emax, max(T_emin, t_e))
      "actual entering condenser water temperature";

    // Chiller capacity fraction
    Unit.Efficiency capFunT;
    Real unit "number of parallel terminals";

    // Chiller energy input ratio
    Unit.Efficiency EIRFunT;

    // Chiller energy input ratio
    Unit.Efficiency EIRFunPLR;

    // Chiller Part Load Ratio (PLR)
    Unit.Efficiency PLR;

    // COP of the chiller
    Unit.Efficiency eta;

    // Total work
    Unit.Power w_tot;


  equation
    // Chiller chiller capacity fraction biquadratic curve
    capFunT = Func.smoothMax(
      x1 = 1E-6,
      x2 = Func.biquadratic(a = per.capFunT, x1 = Tl, x2 = Te),
      deltaX = 1E-7);

    // Chiller energy input ratio biquadratic curve
    EIRFunT = Func.biquadratic(a = per.EIRFunT, x1 = Tl, x2 = Te);

    // Chiller energy input ratio quadratic curve
    //PLR = min(per.PLRMax, max(per.PLRMin, Q/(per.QEva_flow_nominal/unit)/capFunT));
    EIRFunPLR = per.EIRFunPLR[1] + per.EIRFunPLR[2] * PLR +
      per.EIRFunPLR[3] * PLR^2;

    // Calculate COP of the Chiller
    if Q <= deltaQ then
      w_tot = 0.0;
      PLR = 0.0;
      unit = 0.0;
      eta = 0.0;
    else
      eta = per.COP_nominal / EIRFunT / EIRFunPLR * PLR;
      w_tot = Q / eta;
      PLR = 0.5;
      PLR = Q/(-per.QEva_flow_nominal/unit)/capFunT;
    end if;

    annotation (Icon(graphics={
          Rectangle(
            origin={0,1},
            lineColor={64,64,64},
            fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-100.0,-75.0},{100.0,75.0}},
            radius=25.0),
          Text(
            extent={{-6,12},{0,0}},
            lineColor={238,46,47},
            lineThickness=1,
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={217,67,180},
            textString="Chiller COP
Fitting Curve",
            fontSize=10),
          Rectangle(extent={{-100,100},{100,-98}}, lineColor={28,108,200}),
        Text(
          extent={{-160,140},{168,104}},
          lineColor={28,108,200},
          textString="%name")}));
  end ChillerCOPFittingCurve;

  model DesiccantWheel "DW component"

      // Import the package
      package Medium = Modelica.Media.Air.MoistAir;
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput w_set(min=0)
        "supply air humidity ratio"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "process air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealInput t_4(min=0)
        "exhaust air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
      Modelica.Blocks.Interfaces.RealInput rh_4(min=0)
        "exhaust air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput w_c = W_c / 1e3
        "the total work of the desiccant system in kW"
        annotation (Placement(transformation(extent={{100,60},{120,80}}),
          iconTransformation(extent={{100,60},{120,80}})));
      Modelica.Blocks.Interfaces.RealOutput rh_2 = phi_2 * 100
        "supply air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-80},{120,-60}}),
          iconTransformation(extent={{100,-80},{120,-60}})));
      Modelica.Blocks.Interfaces.RealOutput t_2 = Unit.Conversions.to_degC(T_2)
        "supply air temperature in degree C"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
          iconTransformation(extent={{100,-10},{120,10}})));

      // Parameters
      parameter Unit.Pressure P_a = 101325 "air pressure";
      parameter Real theta_o = 1.1
        "dimensionless outlet temperature";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.Efficiency COP_h(min = 0) = 3.0
        "COP of the heater";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in DW to supply air";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the desiccant wheel";
      parameter Unit.Power deltaWc = 1.0 "threshold control on power input";
      parameter Unit.Pressure P_dw "pressure loss across the desiccant wheel";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "process air temperature in K";
      Real phi_1 = rh_1 / 100 "process air relative humidity";
      Real phi_4 = rh_4 / 100 "exhaust air relative humidity";
      Unit.Temperature T_4(displayUnit = "degC") = Unit.Conversions.from_degC(t_4)
        "exhaust air temperature in K";

      // Flow rate of the moist air
      Unit.Density rho_a = Medium.density(State_1) "density of the moist air";
      Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";

      // State 1: Process air
      Real X_1 = Medium.massFraction_pTphi(P_a, T_1, phi_1)
        "mass fraction of water vapor in the process air";
      Medium.ThermodynamicState State_1 = Medium.setState_pTX(P_a, T_1, {X_1})
        "thermodynamic state of process air";
      Real w_1 = X_1 / (1 - X_1) "process air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium.specificEnthalpy(State_1)
        "process air specific enthalpy";

      // State 2: Dehumidified air
      Unit.Temperature T_2(displayUnit = "degC") "supply air temperature";
      Real w_2 "supply air humidity ratio";
      Real X_2 = w_2 / (1 + w_2)
        "mass fraction of water vapor in the supply air";
      Medium.ThermodynamicState State_2 = Medium.setState_pTX(P_a, T_2, {X_2})
        "thermodynamic state of supply air";
      Unit.SpecificEnthalpy h_2 = Medium.specificEnthalpy(State_2)
        "supply air specific enthalpy";
      Real phi_2 = Medium.relativeHumidity(State_2)
        "dehumidified air relative humidity";

      // State 3: Regeneration air
      Unit.Temperature T_3(displayUnit = "degC") "regeneration air temperature";
      Real w_3 = w_4 "regeneration air humidity ratio";
      Real X_3 = X_4 "mass fraction of water vapor in the regeneration air";
      Real phi_3 = phi_is "regeneration air relative humidity";
      Medium.ThermodynamicState State_3 = Medium.setState_pTX(P_a, T_3, {X_3})
        "thermodynamic state of regeneration air";
      Unit.SpecificEnthalpy h_3 = Medium.specificEnthalpy(State_3)
        "regeneration air specific enthalpy";

      // State 4: Exhaust air for regeneration process
      Real X_4 = Medium.massFraction_pTphi(P_a, T_4, phi_4)
        "mass fraction of water vapor in the exhaust air";
      Real w_4 = X_4 / (1 - X_4) "exhaust air humidity ratio";
      Medium.ThermodynamicState State_4 = Medium.setState_pTX(P_a, T_4, {X_4})
        "thermodynamic state of regeneration air";
      Unit.SpecificEnthalpy h_4 = Medium.specificEnthalpy(State_4)
        "regeneration air specific enthalpy";

      // Isenthalpic process
      Unit.SpecificEnthalpy h_is = h_1
        "dehumidified air specific enthalpy for isenthalpic process";
      Real X_is = X_2 "mass fraction of water vapor in the dehumidified air";
      Medium.ThermodynamicState State_is = Medium.setState_phX(P_a, h_is, {X_is})
        "thermodynamic state of dehumidified air for isenthalpic process";
      Unit.Temperature T_is = Medium.temperature(State_is)
        "dehumidified air temperature for isenthalpic process";
      Real phi_is = Medium.relativeHumidity(State_is)
        "dehumidified air relative humidity for isenthalpic process";

      // Heater Power
      Unit.Power W_h "work done by the heater";
      Unit.Power W_p "fan power for the DW";
      Unit.Power W_c "total work done for the DW";


  equation

      // Definition of dimensionless outlet temperature
      // theta_o = (T_2 - T_1) / (T_is - T_1)
      theta_o * (T_is - T_1) = (T_2 - T_1);

      // Calculate the regeneration air temperature
      X_3 = Medium.massFraction_pTphi(P_a, T_3, phi_3);

      // Control on/off the desiccant wheel, only on when latent cooling load is
      // nonzero
      if w_1 - w_set <= deltaw then
        w_2 = w_1;
        W_h = 0.0;
        W_p = 0.0;
        W_c = 0.0;

      else
        w_2 = w_set;
        W_h = max(deltaWc, m_a * r * (h_3 - h_4) / COP_h);
        W_p = P_dw * V_a;
        W_c = max(deltaWc, W_h + W_p);

      end if;


    annotation (Icon(graphics={
          Rectangle(
            extent={{-100,-44},{100,-52}},
            lineColor={238,46,47},
            lineThickness=0.5,
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-100,32},{100,24}},
            lineColor={28,108,200},
            lineThickness=0.5,
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-60,72},{0,-88}},
            lineColor={102,44,145},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{20,-28},{72,-68}},
            lineColor={238,46,47},
            lineThickness=0.5,
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-102,36},{-86,20}},
            lineColor={217,67,180},
            lineThickness=0.5,
            fillColor={217,67,180},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{86,36},{102,20}},
            lineColor={238,46,47},
            lineThickness=0.5,
            fillColor={238,46,47},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-102,-40},{-86,-56}},
            lineColor={28,108,200},
            lineThickness=0.5,
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{86,-40},{102,-56}},
            lineColor={244,125,35},
            lineThickness=0.5,
            fillColor={244,125,35},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{24,-30},{66,-16}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.None,
            textStyle={TextStyle.Bold},
            fontSize=14,
            textString="Heater"),
          Text(
            extent={{-56,72},{-4,86}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.None,
            textStyle={TextStyle.Bold},
            fontSize=14,
            textString="Desiccant Wheel"),
          Text(
            extent={{-164,146},{170,114}},
            lineColor={28,108,200},
            textString="%name")}),           Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Desiccant Wheel</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the desiccant with a heater for regeneration process.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Supply air; State 3: Regeneration air; State 4: Exhaust air for regeneration process</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end DesiccantWheel;


  model MembraneWithCondenser "Membrane Component (w/ condenser)"

      // Import the package
      package Medium1 = Modelica.Media.Air.MoistAir;
      package Medium2 = Modelica.Media.Water.WaterIF97_R2pT; // Steam
      package Medium3 = Modelica.Media.Water.WaterIF97_R1pT; // Liquid
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "process air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
      Modelica.Blocks.Interfaces.RealInput w_set(min=0)
        "supply air humidity ratio"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput q_cond = Q_cond / 1e3
        "condensation cooling load of the membrane system in kW"
        annotation (Placement(transformation(extent={{100,70},{120,90}}),
          iconTransformation(extent={{100,70},{120,90}})));
      Modelica.Blocks.Interfaces.RealOutput w_c = W_tot / 1e3
        "the total work of the membrane system excluding condensation work in kW"
        annotation (Placement(transformation(extent={{100,20},{120,40}}),
          iconTransformation(extent={{100,20},{120,40}})));
      Modelica.Blocks.Interfaces.RealOutput t_s = Unit.Conversions.to_degC(T_2)
        "supply air temperature in degree C"
        annotation (Placement(transformation(extent={{100,-40},{120,-20}}),
          iconTransformation(extent={{100,-40},{120,-20}})));
      Modelica.Blocks.Interfaces.RealOutput rh_s = phi_2 * 100
        "supply air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-90},{120,-70}}),
          iconTransformation(extent={{100,-90},{120,-70}})));

      // Parameters
      parameter Unit.Pressure P_a = 101325 "air pressure";
      parameter Unit.Efficiency eta_is(max = 1.0) = 0.9
        "isentropic efficiency of the vacuum pump";
      parameter Unit.Pressure P_stage = 2000
        "pressure at the exit of stage pump";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Real deltaw_set(min = 0) = 1e-5
        "humidity threshold control on/off for the membrane system";
      parameter Real deltaP_stage(min = 0) = 1
        "pressure threshold control on the vacuum pump";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on power input";
      parameter Unit.Area A_m = 1.0
        "Membrane area per channel in m2";
      parameter Real n_m = 100.0 "number of membrane stacks";
      parameter Unit.Length h = 0.01 "Channel height in m";
      parameter Real p_w = 6.2e-7 "water permeance of the membrane in mol/m2 s Pa";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "process air temperature";
      Real phi_1 = rh_1 / 100 "process air relative humidity";

      // Flow rate of the moist air
      Unit.Density rho_a = Medium1.density(State_1) "density of the moist air";
      Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";
      Unit.Velocity v_a "velocity of air in m/s";
      Unit.DynamicViscosity mu_a = Medium1.dynamicViscosity(State_1)
        "viscosity of air";
      Unit.Length L = sqrt(A_m) "channel length in m";

      // Flow rate of water vapor
      Unit.Density rho_w = Medium3.density(State_5) "density of liquid water";
      Unit.MassFlowRate m_w "mass flow rate of water";
      Unit.VolumeFlowRate V_w "volumetric flow rate of liquid water";
      Real mw_w = 0.018 "molecular weight of water in kg/mol";
      Real deltaw "humidity removal target for the membrane dehumidifier";

      // State 1: Process air
      Real X_1 = Medium1.massFraction_pTphi(P_a, T_1, phi_1)
        "mass fraction of water vapor in the process air";
      Medium1.ThermodynamicState State_1 = Medium1.setState_pTX(P_a, T_1, {X_1})
        "thermodynamic state of process air";
      Real w_1 = X_1 / (1 - X_1) "process air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium1.specificEnthalpy(State_1)
        "process air specific enthalpy";
      Unit.Pressure P_1v = Bldg.pW_X(X_1) "water vapor pressure at State 1";

      // State 2: Supply air
      Unit.Temperature T_2(displayUnit = "degC") = T_1 "supply air temperature";
      Real w_2 = w_set "supply air humidity ratio";
      Real X_2 = w_2 / (1 + w_2)
        "mass fraction of water vapor in the supply air";
      Medium1.ThermodynamicState State_2 = Medium1.setState_pTX(P_a, T_2, {X_2})
        "thermodynamic state of supply air";
      Unit.SpecificEnthalpy h_2 = Medium1.specificEnthalpy(State_2)
        "supply air specific enthalpy";
      Real phi_2 = Medium1.relativeHumidity(State_2)
        "dehumidified air relative humidity";
      Unit.Pressure P_2v = Bldg.pW_X(X_2) "water vapor pressure at State 2";

      // State 3: Water vapor at the sweep side of the membrane
      Real X_3 = 1.0 "mass fraction of water vapor in the sweep-side steam";
      Unit.Temperature T_3(displayUnit = "degC") = T_1
        "sweep-side water vapor temperature";
      Unit.Pressure P_3 = P_sweep "sweep-side water vapor pressure";
      Medium2.ThermodynamicState State_3 = Medium2.setState_pTX(P_3, T_3, {X_3});
      Unit.SpecificEnthalpy h_3 = Medium2.specificEnthalpy(State_3)
        "sweep-side water vapor specific enthalpy";
      Unit.SpecificEntropy s_3 = Medium2.specificEntropy(State_3)
        "sweep-side water vapor specific entropy";

      // State 4: Pressurized water vapor at the exit of vacuum pump
      Real X_4 = 1.0 "mass fraction of pressurized water vapor";
      Unit.Pressure P_4 = if P_stage >= P_sweep then max(P_1v, P_stage)
        else max(P_1v, P_sweep + deltaP_stage) "pressurized water vapor pressure";
      // Isentropic process from State 3, use the same thermodynamic record
      //   of specific entropy
      Unit.SpecificEnthalpy h_is = Medium2.isentropicEnthalpy(P_4, State_3)
        "isentropic specific enthalpy";
      Unit.SpecificEnthalpy h_4 "pressurized water vapor specific enthalpy";
      Medium2.ThermodynamicState State_4 = Medium2.setState_phX(P_4, h_4, {X_4});
      Unit.Temperature T_4(displayUnit = "degC") = Medium2.temperature(State_4)
        "pressurized water vapor temperature";

      // State 5: Pressurized liquid water (Condensate)
      Real X_5 = 1.0 "mass fraction of liquid water in the condensate";
      Unit.Pressure P_5 = P_4 "condensate pressure";
      // Condensate temperature = Saturation temperature of water vapor
      Unit.Temperature T_5(displayUnit = "degC")=
        Medium2.saturationTemperature(P_4) "pressurized liquid temperature";
      Medium3.ThermodynamicState State_5 = Medium3.setState_pTX(P_5, T_5, {X_5});
      Unit.SpecificEnthalpy h_5 = Medium3.specificEnthalpy(State_5)
        "liquid water specific enthalpy";

      // Calculate the COP of the membrane system
      Unit.Power W_stage "work done by the stage vacuum pump";
      Unit.Power W_pump "work done by the liquid compressor";
      Unit.Power Q_cond "heat to be removed by the condenser";
      Unit.Power W_tot "the total work of the membrane sytem";

      // Pressure loss across the membrane
      Unit.Pressure P_sweep(min = 0)
        "pump pressure at the sweep side of membrane in Pa";
      Unit.Pressure P_avg(min = 0) "average feed pressure along the membrane in Pa";
      Unit.Power W_c "work done by the pump";
      Unit.Pressure deltaP "pressure drop along the pipe";
      Real f "friction factor";
      Real pi = Modelica.Constants.pi;
      Real Re "Reynolds number of pipe flow";

  equation
      // Average water vapor partial pressure across the membrane
      P_avg = (P_1v - P_2v) / Modelica.Math.log(P_1v / P_2v);

      // Volumetric flow rate of liquid water
      deltaw = w_1 - w_set;
      m_w = m_a * deltaw;
      V_w = m_w / rho_w;
      m_w = p_w * mw_w * A_m * n_m * (P_avg - P_sweep);

      // Pressure drop across the membrane
      deltaP = f * 1/2 * rho_a * v_a^2 * L / h;
      m_a = rho_a * v_a * (n_m * L * h);
      Re = rho_a * v_a * 2 * h / mu_a;
      if (Re <= 0) then
        f = 0;
      else
        if (Re < 2.3e3) then
          f = 64 / Re;
        else
          f = (1.82 * Modelica.Math.log10(Re) - 1.64) ^ (-2);
        end if;
      end if;
      W_c = deltaP * V_a;

      // Work done by stage vacuum pump
      eta_is = (h_is - h_3) / (h_4 - h_3);
      W_stage = m_w * (h_4 - h_3);

      // Work done by condenser
      Q_cond = m_w * (h_4 - h_5);

      // Work done by liquid compressor
      W_pump = (P_a - P_4) * V_w;

      // COP of the membrane system
      if deltaw <= deltaw_set then
        W_tot = 0.0;
      else
        W_tot = max(deltaWc, W_stage + W_pump + W_c);
      end if;

    annotation (Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-80,60},{40,40}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-80,40},{-20,20},{40,40},{-80,40}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-10,10},{10,-10}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,127,0},
            fillPattern=FillPattern.Solid,
            origin={-22,-6},
            rotation=-90),
          Polygon(
            points={{-10,4},{10,4},{10,-4},{0,-4},{-10,4}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,127,0},
            fillPattern=FillPattern.Solid,
            origin={-16,-16},
            rotation=-90),
          Line(
            points={{-20,20},{-20,4}},
            color={0,0,0},
            thickness=0.5),
          Rectangle(
            extent={{-40,-40},{0,-60}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-40,-50},{0,-50}},
            color={0,0,0},
            thickness=0.5,
            pattern=LinePattern.Dash),
          Line(
            points={{-16,-26},{-16,-40}},
            color={0,0,0},
            thickness=0.5),
          Ellipse(
            extent={{20,-40},{40,-60}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,140,72},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{30,-40},{50,-40},{50,-48},{40,-48},{30,-40}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,140,72},
            fillPattern=FillPattern.Solid),
          Line(
            points={{0,-56},{22,-56}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{50,-44},{84,-44}},
            color={0,0,0},
            thickness=0.5,
            arrow={Arrow.None,Arrow.Filled}),
          Text(
            extent={{-78,-64},{38,-68}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=20,
            fontName="serif",
            textString="Condenser"),
          Text(
            extent={{-78,68},{38,64}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=20,
            fontName="serif",
            textString="Membrane"),
        Text(
          extent={{-160,140},{168,108}},
          lineColor={28,108,200},
          textString="%name"),
          Text(
            extent={{-96,48},{-78,46}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="1"),
          Text(
            extent={{38,48},{56,46}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="2"),
          Text(
            extent={{-22,12},{-4,10}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="3"),
          Text(
            extent={{-20,-32},{-2,-34}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="4"),
          Text(
            extent={{2,-50},{20,-52}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="5")}),        Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Membrane Dehumidifier with Condenser Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the membrane dehumidifier with a condenser recycling excess water vapor in the air.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Supply air; State 3: Water vapor at the sweep side of the membrane; State 4: Pressurized water vapor for condensation; State 5: Pressurized liquid water for reuse</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end MembraneWithCondenser;





  model MembraneWithSweepGas

      // Import the package
      package Medium1 = Modelica.Media.Air.MoistAir;
      package Medium2 = Modelica.Media.Water.WaterIF97_R2pT; // Steam
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "process air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealInput w_set(min=0)
        "supply air humidity ratio"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_3(min=0)
        "exhaust air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
      Modelica.Blocks.Interfaces.RealInput rh_3(min=0)
        "exhaust air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput w_c = W_c / 1e3
        "the total work of the membrane system in kW"
        annotation (Placement(transformation(extent={{100,70},{120,90}}),
          iconTransformation(extent={{100,70},{120,90}})));
      Modelica.Blocks.Interfaces.RealOutput rh_s = phi_2 * 100
        "process air relative humidity IN %"
        annotation (Placement(transformation(extent={{100,-90},{120,-70}}),
          iconTransformation(extent={{100,-90},{120,-70}})));
      Modelica.Blocks.Interfaces.RealOutput t_s = Unit.Conversions.to_degC(T_2)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
          iconTransformation(extent={{100,-10},{120,10}})));


      // Parameters
      parameter Unit.Pressure P_a = 101325 "air pressure";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in DW to supply air";
      parameter Real deltaw_set(min = 0) = 1e-5
        "humidity threshold control on/off for the membrane system";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on power input";
      parameter Unit.Area A_m = 1.0
        "Membrane area per channel in m2";
      parameter Real n_m = 100.0 "number of membrane stacks";
      parameter Unit.Length h = 0.01 "Channel height in m";
      parameter Real p_w = 6.2e-7 "water permeance of the membrane in mol/m2 s Pa";
      parameter Unit.Efficiency COP_h = 4.8 "COP of heat pump";
      parameter Unit.Efficiency eta_is(max = 1.0) = 0.9
        "isentropic efficiency of the vacuum pump";


      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "process air temperature";
      Real phi_1 = rh_1 / 100 "process air relative humidity";
      Unit.Temperature T_3(displayUnit = "degC") = Unit.Conversions.from_degC(t_3)
        "exhaust air temperature";
      Real phi_3 = rh_3 / 100 "exhaust air relative humidity";

      // Flow rate of the moist air
      Unit.Density rho_a = Medium1.density(State_1) "density of the moist air";
      Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";
      Unit.Velocity v_a "velocity of air in m/s";
      Unit.DynamicViscosity mu_a = Medium1.dynamicViscosity(State_1)
        "viscosity of air";
      Unit.Length L = sqrt(A_m) "channel length in m";

      // Flow rate of water vapor
      Unit.MassFlowRate m_w "mass flow rate of water";
      Real mw_w = 0.018 "molecular weight of water in kg/mol";
      Real deltaw "humidity removal target for the membrane dehumidifier";

      // State 1: Process air
      Real X_1 = Medium1.massFraction_pTphi(P_a, T_1, phi_1)
        "mass fraction of water vapor in the process air";
      Medium1.ThermodynamicState State_1 = Medium1.setState_pTX(P_a, T_1, {X_1})
        "thermodynamic state of process air";
      Real w_1 = X_1 / (1 - X_1) "process air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium1.specificEnthalpy(State_1)
        "process air specific enthalpy";
      Unit.Pressure P_1v = Bldg.pW_X(X_1) "water vapor pressure at State 1";

      // State 2: Supply air
      Unit.Temperature T_2(displayUnit = "degC") "supply air temperature";
      Real w_2 = w_set "supply air humidity ratio";
      Real X_2 = w_2 / (1 + w_2)
        "mass fraction of water vapor in the supply air";
      Unit.Pressure P_2v = Bldg.pW_X(X_2) "water vapor pressure at State 2";
      Medium1.ThermodynamicState State_2 = Medium1.setState_pTX(P_a, T_3, {X_2})
        "thermodynamic state of supply air";
      Unit.SpecificEnthalpy h_2 = Medium1.specificEnthalpy(State_2)
        "supply air specific enthalpy";
      Real phi_2 = Medium1.relativeHumidity(State_2)
        "supply air relative humidity";

      // State 3: Exhaust air
      Real X_3 = Medium1.massFraction_pTphi(P_a, T_3, phi_3)
        "mass fraction of water vapor in the exhaust air";
      Medium1.ThermodynamicState State_3 = Medium1.setState_pTX(P_a, T_3, {X_3})
        "thermodynamic state of exhaust air";
      Real w_3 = X_3 / (1 - X_3) "exhaust air humidity ratio";
      Unit.SpecificEnthalpy h_3 = Medium1.specificEnthalpy(State_3)
        "exhaust air specific enthalpy";

      // State 4: Sweeping air
      Real X_4 = X_3 "mass fraction of water vapor in the sweeping air";
      Real w_4 = X_4 / (1 - X_4) "sweeping air humidity ratio";
      Real P_4(max = 1e5) "exhaust air pressure after throttling";
      Unit.Pressure P_4v "water vapor pressure at State 4";
      Unit.Temperature T_4(displayUnit = "degC") = max(T_4min1, max(T_4min2, T_4min3))
        "sweeping air temperature";
      Medium1.ThermodynamicState State_4 = Medium1.setState_pTX(P_4, T_4, {X_4})
        "thermodynamic state of sweeping air";
      Unit.SpecificEnthalpy h_4 = Medium1.specificEnthalpy(State_4)
        "sweeping air specific enthalpy";

      // Sweeping air after throttling
      // Using throtting valve to reduce the temperature and pressure
      // Note that the throtting process is isenthalpic and adiabatic, and thus
      // it has no work input. (h_4min1 = h_3)
      Unit.Temperature T_4min1(displayUnit = "degC") = Medium1.T_phX(P_4, h_3, {X_3})
        "sweeping air temperature after throttling";
      Medium1.ThermodynamicState State_4min1 = Medium1.setState_pTX(P_4, T_4, {X_4})
        "thermodynamic state of sweeping air after throttling";

      // Sweeping air minimum temperature determined by condensation at State 4
      Unit.Temperature T_4min2 = Bldg.TDewPoi_pW(P_4v)
        "sweeping air minimum temperature due to condensation";

      // Sweeping air minimum temperature determined by condensation at State 5
      Unit.SpecificEnthalpy h_4min3 "sweeping air specific enthalpy";
      Medium1.ThermodynamicState State_4min3 = Medium1.setState_phX(P_4, h_4min3, {X_4})
        "thermodynamic state of sweeping air with minimum temperature";
      Unit.Temperature T_4min3(displayUnit = "degC") = Medium1.temperature(State_4min3)
        "sweeping air minimum temperature to prevent condensation in the exhaust";

      // State 5min: Exhaust sweeping air with minimum temperature
      Real w_5 = w_4 + deltaw "exhaust sweeping air humidity ratio";
      Real X_5 = w_5 / (1 + w_5)
        "mass fraction of water vapor in the exhaust sweeping air";
      Unit.Pressure P_5v = Bldg.pW_X(X_5, p = P_4) "water vapor pressure at State 5";
      Unit.Temperature T_5min(displayUnit = "degC") = Bldg.TDewPoi_pW(P_5v)
        "dewpoint temperature of exhaust sweeping air";
      Medium1.ThermodynamicState State_5min = Medium1.setState_pTX(P_4, T_5min, {X_5})
        "minimum thermodynamic state of exhaust sweeping air";
      Unit.SpecificEnthalpy h_5min = Medium1.specificEnthalpy(State_5min)
        "exhaust sweeping air minimum specific enthalpy";

      // State 5: Exhaust sweeping air (P_5 = P_4)
      Medium1.ThermodynamicState State_5 = Medium1.setState_phX(P_4, h_5, {X_5})
        "thermodynamic state of exhaust sweeping air";
      Unit.Temperature T_5(displayUnit = "degC") = Medium1.temperature(State_5)
        "exhaust sweeping air temperature";
      Unit.SpecificEnthalpy h_5 "exhaust sweeping air specific enthalpy";

      // State 6: Pressurized exhaust air
      Unit.Pressure P_6 = 40000 "pressurized exhaust air pressure";
      Real X_6 = X_5 "mass fraction of water vapor in the pressurized exhaust air";
      // Isentropic process from State 5, use the same thermodynamic record
      //   of specific entropy
      Unit.SpecificEnthalpy h_is = Medium1.isentropicEnthalpyApproximation(P_6, State_5)
        "isentropic specific enthalpy";
      Unit.SpecificEnthalpy h_6 "pressurized exhaust air specific enthalpy";
      Medium1.ThermodynamicState State_6 = Medium1.setState_phX(P_6, h_6, {X_5});
      Unit.Temperature T_6(displayUnit = "degC") = Medium1.temperature(State_6)
        "pressurized moist air temperature";
      Real phi_6 = Medium1.relativeHumidity(State_6);

      // Pressure loss across the membrane
      Unit.Power W_c "work done by the pump";
      Unit.Pressure deltaP "pressure drop along the pipe";
      Real f "friction factor";
      Real pi = Modelica.Constants.pi;
      Real Re "Reynolds number of pipe flow";
      Unit.Efficiency epsilon
        "effectiveness of sensible heat recovery across the membrane";
      Real NTU = A_m * n_m * 10 / m_a / 1000
        "NTU of heat exchange in the membrane, assuming U = 10 W/m2 K and Cp = 1000";


  equation
      // Volumetric flow rate of liquid water
      deltaw = w_1 - w_set;
      m_w = m_a * deltaw;
      m_w = p_w * mw_w * A_m * n_m * (P_2v - P_4v);

      // Find the expansion pressure
      X_4 = Bldg.X_pW(P_4v, p = P_4);

      // Heat balance of the membrane
      epsilon = (T_1 - T_2) / (T_1 - T_4);
      epsilon = NTU / (1 + NTU);
      h_1 + r * h_4min3 = h_2 + r * h_5min;
      h_1 + r * h_4 = h_2 + r * h_5;

      // Pressure drop across the membrane
      deltaP = f * 1/2 * rho_a * v_a^2 * L / h;
      m_a = rho_a * v_a * (n_m * L * h);
      Re = rho_a * v_a * 2 * h / mu_a;
      if (Re <= 0) then
        f = 0;
      else
        if (Re < 2.3e3) then
          f = 64 / Re;
        else
          f = (1.82 * Modelica.Math.log10(Re) - 1.64) ^ (-2);
        end if;
      end if;

      // Work done by vacuum pump
      eta_is = (h_is - h_5) / (h_6 - h_5);

      // COP of the membrane system
      if deltaw <= deltaw_set then
        W_c = 0.0;
      else
        W_c = max(deltaWc,
                  deltaP * V_a + max(0, m_a * (h_4 - h_3) / COP_h) + m_a * (h_6 - h_5));
      end if;

    annotation (Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-80,20},{40,0}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-78,32},{38,28}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=20,
            fontName="serif",
            textString="Membrane"),
        Text(
          extent={{-160,140},{168,108}},
          lineColor={28,108,200},
          textString="%name"),
          Text(
            extent={{-96,10},{-78,8}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="1"),
          Text(
            extent={{38,10},{56,8}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="2"),
          Rectangle(
            extent={{-80,0},{40,-20}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{38,-6},{56,-8}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="4"),
          Text(
            extent={{-96,-10},{-78,-12}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="5"),
          Line(
            points={{56,-8},{76,-18},{76,-8},{56,-18},{56,-16},{56,-8}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{40,-12},{56,-12}},
            color={0,0,0},
            thickness=0.5,
            arrow={Arrow.Filled,Arrow.None}),
          Text(
            extent={{76,-6},{94,-8}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="3'")}),       Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Membrane Dehumidifier with Condenser Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the membrane dehumidifier with a condenser recycling excess water vapor in the air.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Supply air; State 3: Water vapor at the sweep side of the membrane; State 4: Pressurized water vapor for condensation; State 5: Pressurized liquid water for reuse</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end MembraneWithSweepGas;

  model MembraneWithMembrane

      // Import the package
      package Medium1 = Modelica.Media.Air.MoistAir;
      package Medium2 = Modelica.Media.Water.WaterIF97_R2pT; // Steam
      package Unit = Modelica.SIunits;
      package Bldg = Buildings.Utilities.Psychrometrics.Functions;

      // Inputs
      Modelica.Blocks.Interfaces.RealInput t_1(min=0)
        "process air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
      Modelica.Blocks.Interfaces.RealInput rh_1(min=0)
        "process air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Blocks.Interfaces.RealInput w_set(min=0)
        "supply air humidity ratio"
        annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
      Modelica.Blocks.Interfaces.RealInput t_5(min=0)
        "exhaust air temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
      Modelica.Blocks.Interfaces.RealInput rh_5(min=0)
        "exhaust air relative humidity in %"
        annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));

      // Outputs
      Modelica.Blocks.Interfaces.RealOutput w_c = W_tot / 1e3
        "the total work of the membrane system excluding condensation work in kW"
        annotation (Placement(transformation(extent={{100,70},{120,90}}),
          iconTransformation(extent={{100,70},{120,90}})));
      Modelica.Blocks.Interfaces.RealOutput t_s = Unit.Conversions.to_degC(T_2)
        "supply air temperature in degree C"
        annotation (Placement(transformation(extent={{100,-10},{120,10}}),
          iconTransformation(extent={{100,-10},{120,10}})));
      Modelica.Blocks.Interfaces.RealOutput rh_s = phi_2 * 100
        "supply air relative humidity in %"
        annotation (Placement(transformation(extent={{100,-90},{120,-70}}),
          iconTransformation(extent={{100,-90},{120,-70}})));

      // Parameters
      parameter Unit.Pressure P_a = 101325 "air pressure";
      parameter Unit.Efficiency eta_is(max = 1.0) = 0.9
        "isentropic efficiency of the vacuum pump";
      //parameter Unit.Pressure P_stage = 2000
      //  "pressure at the exit of stage pump";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Real deltaw_set(min = 0) = 1e-5
        "humidity threshold control on/off for the membrane system";
      parameter Real deltaP_stage(min = 0) = 1
        "pressure threshold control on the vacuum pump";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on power input";
      parameter Unit.Area A_m = 1.0
        "Membrane area per channel in m2";
      parameter Real n_m = 100.0 "number of membrane stacks";
      parameter Unit.Length h = 0.01 "Channel height in m";
      parameter Real p_w = 6.2e-7 "water permeance of the membrane in mol/m2 s Pa";
      parameter Unit.Efficiency COP_h = 4.8 "COP of heat pump";

      // Conversion of inputs to physical parameters
      Unit.Temperature T_1(displayUnit = "degC") = Unit.Conversions.from_degC(t_1)
        "process air temperature";
      Real phi_1 = rh_1 / 100 "process air relative humidity";
      Unit.Temperature T_5(displayUnit = "degC") = Unit.Conversions.from_degC(t_5)
        "exhaust air temperature";
      Real phi_5 = rh_5 / 100 "exhaust air relative humidity";

      // Flow rate of the moist air
      Unit.Density rho_a = Medium1.density(State_1) "density of the moist air";
      Unit.MassFlowRate m_a = rho_a * V_a "mass flow rate of the moist air";
      Unit.Velocity v_a "velocity of air in m/s";
      Unit.DynamicViscosity mu_a = Medium1.dynamicViscosity(State_1)
        "viscosity of air";
      Unit.Length L = sqrt(A_m) "channel length in m";

      // Flow rate of water vapor
      Unit.MassFlowRate m_w "mass flow rate of water";
      Real mw_w = 0.018 "molecular weight of water in kg/mol";
      Real deltaw "humidity removal target for the membrane dehumidifier";

      // State 1: Process air
      Real X_1 = Medium1.massFraction_pTphi(P_a, T_1, phi_1)
        "mass fraction of water vapor in the process air";
      Medium1.ThermodynamicState State_1 = Medium1.setState_pTX(P_a, T_1, {X_1})
        "thermodynamic state of process air";
      Real w_1 = X_1 / (1 - X_1) "process air humidity ratio";
      Unit.SpecificEnthalpy h_1 = Medium1.specificEnthalpy(State_1)
        "process air specific enthalpy";
      Unit.Pressure P_1v = Bldg.pW_X(X_1) "water vapor pressure at State 1";

      // State 2: Supply air
      Unit.Temperature T_2(displayUnit = "degC") = T_1 "supply air temperature";
      Real w_2 = w_set "supply air humidity ratio";
      Real X_2 = w_2 / (1 + w_2)
        "mass fraction of water vapor in the supply air";
      Medium1.ThermodynamicState State_2 = Medium1.setState_pTX(P_a, T_2, {X_2})
        "thermodynamic state of supply air";
      Unit.SpecificEnthalpy h_2 = Medium1.specificEnthalpy(State_2)
        "supply air specific enthalpy";
      Real phi_2 = Medium1.relativeHumidity(State_2)
        "dehumidified air relative humidity";
      Unit.Pressure P_2v = Bldg.pW_X(X_2) "water vapor pressure at State 2";

      // State 3: Water vapor at the sweep side of the membrane
      Real X_3 = 1.0 "mass fraction of water vapor in the sweep-side steam";
      Unit.Temperature T_3(displayUnit = "degC") = T_1
        "sweep-side water vapor temperature";
      Unit.Pressure P_3 = P_sweep "sweep-side water vapor pressure";
      Medium2.ThermodynamicState State_3 = Medium2.setState_pTX(P_3, T_3, {X_3});
      Unit.SpecificEnthalpy h_3 = Medium2.specificEnthalpy(State_3)
        "sweep-side water vapor specific enthalpy";
      Unit.SpecificEntropy s_3 = Medium2.specificEntropy(State_3)
        "sweep-side water vapor specific entropy";

      // State 4: Pressurized water vapor at the exit of vacuum pump
      Real X_4 = 1.0 "mass fraction of pressurized water vapor";
      Unit.Pressure P_4 = if P_stage >= P_sweep then P_stage
        else P_sweep + deltaP_stage "pressurized water vapor pressure";
      // Isentropic process from State 3, use the same thermodynamic record
      //   of specific entropy
      Unit.SpecificEnthalpy h_is = Medium2.isentropicEnthalpy(P_4, State_3)
        "isentropic specific enthalpy";
      Unit.SpecificEnthalpy h_4 "pressurized water vapor specific enthalpy";
      Medium2.ThermodynamicState State_4 = Medium2.setState_phX(P_4, h_4, {X_4});
      Unit.Temperature T_4(displayUnit = "degC") = Medium2.temperature(State_4)
        "pressurized water vapor temperature";

      // State 5: Exhaust air
      Real X_5 = Medium1.massFraction_pTphi(P_a, T_5, phi_5)
        "mass fraction of water vapor in the exhaust air";
      Medium1.ThermodynamicState State_5 = Medium1.setState_pTX(P_a, T_5, {X_5})
        "thermodynamic state of exhaust air";
      Real w_5 = X_5 / (1 - X_5) "exhaust air humidity ratio";
      Unit.SpecificEnthalpy h_5 = Medium1.specificEnthalpy(State_5)
        "exhaust air specific enthalpy";
      Unit.Pressure P_5v = Bldg.pW_X(X_5) "water vapor pressure at State 5";

      // State 6: Humidified exhaust air
      Unit.SpecificEnthalpy h_6 "humidified exhaust air specific enthalpy";
      Real w_6 = w_5 + deltaw "humidified exhaust air humidity ratio";
      Real X_6 = w_6 / (1 + w_6)
        "mass fraction of water vapor in the humidified exhaust air";
      Medium1.ThermodynamicState State_6 = Medium1.setState_phX(P_a, h_6, {X_6})
        "thermodynamic state of humidified exhaust air";
      Unit.Temperature T_6(displayUnit = "degC") = Medium1.temperature(State_6)
        "humidified exhaust air temperature";
      Unit.Pressure P_6v = Bldg.pW_X(X_6) "water vapor pressure at State 6";

      // State 5': Heated exhaust air to avoid condensation at State 6
      Unit.Temperature T_6min(displayUnit = "degC") = Bldg.TDewPoi_pW(P_6v)
        "dewpoint temperature of humidified exhaust air";
      Medium1.ThermodynamicState State_6min = Medium1.setState_pTX(P_a, T_6min, {X_6})
        "minimum thermodynamic state of humidified exhaust air";
      Unit.SpecificEnthalpy h_6min = Medium1.specificEnthalpy(State_6min)
        "humidified exhaust air minimum specific enthalpy";
      Unit.SpecificEnthalpy h_5min "exhaust air minimum specific enthalpy";
      Medium1.ThermodynamicState State_5min = Medium1.setState_phX(P_a, h_5min, {X_5})
        "minimum thermodynamic state of exhaust air";
      Unit.Temperature T_5min(displayUnit = "degC") = Medium1.temperature(State_5min)
        "humidified exhaust air minimum temperature";

      // Calculate the COP of the membrane system
      Unit.Power W_stage "work done by the stage vacuum pump";
      Unit.Power W_h "work done by the heat pump";
      Unit.Power W_tot "the total work of the membrane sytem";

      // Pressure loss across the membrane
      Unit.Pressure P_sweep(min = 0)
        "pump pressure at the sweep side of membrane in Pa";
      Unit.Pressure P_stage(min = 0) "pump pressure at the exit in Pa";
      Unit.Pressure P_avg_p(min = 0)
        "average process pressure along the membrane in Pa";
      Unit.Pressure P_avg_e(min = 0)
        "average exhaust air pressure along the membrane in Pa";
      Unit.Power W_c "work done by the pump";
      Unit.Pressure deltaP "pressure drop along the pipe";
      Real f "friction factor";
      Real pi = Modelica.Constants.pi;
      Real Re "Reynolds number of pipe flow";

  equation
      // Average water vapor partial pressure across the membrane
      P_avg_p = (P_1v - P_2v) / Modelica.Math.log(P_1v / P_2v);
      P_avg_e = (P_6v - P_5v) / Modelica.Math.log(P_6v / P_5v);

      // Volumetric flow rate of liquid water
      deltaw = w_1 - w_set;
      m_w = m_a * deltaw;
      m_w = p_w * mw_w * A_m * n_m * (P_avg_p - P_sweep);

      // Mass balance of water content in the system
      P_avg_p - P_sweep = P_stage - P_avg_e;

      // Heat transfer from hot water vapor to the exhaust air
      m_w * h_4 + m_a * h_5 = m_a * h_6;
      m_w * h_4 + m_a * h_5min = m_a * h_6min;

      // Pressure drop across the membrane
      deltaP = f * 1/2 * rho_a * v_a^2 * L / h;
      m_a = rho_a * v_a * (n_m * L * h);
      Re = rho_a * v_a * 2 * h / mu_a;
      if (Re <= 0) then
        f = 0;
      else
        if (Re < 2.3e3) then
          f = 64 / Re;
        else
          f = (1.82 * Modelica.Math.log10(Re) - 1.64) ^ (-2);
        end if;
      end if;
      W_c = deltaP * V_a;

      // Work done by stage vacuum pump
      eta_is = (h_is - h_3) / (h_4 - h_3);
      W_stage = m_w * (h_4 - h_3);

      // Work done by the heat pump
      W_h = max(0, m_a * (h_5min - h_5) / COP_h);

      // COP of the membrane system
      if deltaw <= deltaw_set then
        W_tot = 0.0;
      else
        W_tot = max(deltaWc, W_stage + W_h + 2 * W_c);
      end if;

    annotation (Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
          Rectangle(
            extent={{-80,70},{40,50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-80,50},{-20,30},{40,50},{-80,50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-10,10},{10,-10}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,127,0},
            fillPattern=FillPattern.Solid,
            origin={-22,4},
            rotation=-90),
          Polygon(
            points={{-10,4},{10,4},{10,-4},{0,-4},{-10,4}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={0,127,0},
            fillPattern=FillPattern.Solid,
            origin={-16,-6},
            rotation=-90),
          Line(
            points={{-20,30},{-20,14}},
            color={0,0,0},
            thickness=0.5),
          Line(
            points={{-16,-16},{-16,-30}},
            color={0,0,0},
            thickness=0.5),
          Text(
            extent={{-78,78},{38,74}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=20,
            fontName="serif",
            textString="Membrane"),
        Text(
          extent={{-160,140},{168,108}},
          lineColor={28,108,200},
          textString="%name"),
          Text(
            extent={{-96,58},{-78,56}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="1"),
          Text(
            extent={{38,58},{56,56}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="2"),
          Text(
            extent={{-22,22},{-4,20}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="3"),
          Text(
            extent={{-20,-22},{-2,-24}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="4"),
          Text(
            extent={{42,-58},{60,-60}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="5"),
          Rectangle(
            extent={{-76,-70},{44,-50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-76,-50},{-16,-30},{44,-50},{-76,-50}},
            lineColor={0,0,0},
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-74,-80},{42,-76}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=20,
            fontName="serif",
            textString="Membrane"),
          Text(
            extent={{-92,-60},{-74,-62}},
            lineColor={0,0,0},
            pattern=LinePattern.Dash,
            lineThickness=0.5,
            fillColor={102,44,145},
            fillPattern=FillPattern.Solid,
            fontSize=10,
            fontName="serif",
            textString="6")}),        Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Membrane Dehumidifier with Condenser Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of the membrane dehumidifier with a condenser recycling excess water vapor in the air.</span></p>
<p><br><span style=\"font-family: Arial;\">State 1: Hot air inlet; State 2: Supply air; State 3: Water vapor at the sweep side of the membrane; State 4: Pressurized water vapor for condensation; State 5: Pressurized liquid water for reuse</span></p>
<p><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end MembraneWithMembrane;

  model MinWork

    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";

      // Control Parameters
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Efficiency epsilon_s = 0.80 "sensible heat recovery ratio";
      parameter Unit.Efficiency epsilon_l = 0.77 "latent heat recovery ratio";
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,64},{120,84}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-40,-20},{-20,0}})));

    HVACThermo.MinWork4DOASCooling minWork4DOASCooling(P_1 = P, P_2 = P, V_a = V_a,
      epsilon_s = epsilon_s, epsilon_l = epsilon_l, deltaT = deltaT, deltaw = deltaw,
      deltaP = deltaP)
      annotation (Placement(transformation(extent={{20,-20},{40,0}})));
    Modelica.Blocks.Interfaces.RealOutput w_min_sa
      "Minimum work for saturated air in kW"
      annotation (Placement(transformation(extent={{100,-18},{120,2}})));
    Modelica.Blocks.Interfaces.RealOutput w_min_l
      "minimum work for liquid water in kW"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.RealOutput w_min_erw
      "minimum work for DOAS cooling system in kW if ERW is implemented"
      annotation (Placement(transformation(extent={{100,-84},{120,-64}})));
    Modelica.Blocks.Interfaces.RealOutput w_min_hw
      "minimum work for DOAS cooling system in kW if HW is implemented"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat)
      annotation (Line(points={{-120,80},{-42,80},{-42,-1}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-88,
            40},{-88,-5.4},{-42,-5.4}},
                                    color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-100,0},
            {-100,-10},{-42,-10}},
                                color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3)
      annotation (Line(points={{-120,-80},{-42,-80},{-42,-19}},color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{
            -92,-40},{-92,-14.6},{-42,-14.6}},
                                    color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-19,-4},{-11.5,
            -4},{-11.5,74},{110,74}},
                                  color={0,0,127}));
    connect(dOASCoolingLoad.w_set, minWork4DOASCooling.w_set) annotation (Line(
          points={{-19,-10},{-6,-10},{-6,-2},{18,-2}}, color={0,0,127}));
    connect(t_o, minWork4DOASCooling.t_1) annotation (Line(points={{-120,40},{8,40},
            {8,-6},{18,-6}}, color={0,0,127}));
    connect(rh_i, minWork4DOASCooling.rh_4) annotation (Line(points={{-120,-80},{8,
            -80},{8,-18},{18,-18}}, color={0,0,127}));

    connect(rh_o, minWork4DOASCooling.rh_1) annotation (Line(points={{-120,-40},{0,
            -40},{0,-10},{18,-10}}, color={0,0,127}));
    connect(t_i, minWork4DOASCooling.t_4) annotation (Line(points={{-120,0},{-100,
            0},{-100,-10},{-60,-10},{-60,-30},{4,-30},{4,-14},{18,-14}}, color={0,
            0,127}));
    connect(minWork4DOASCooling.w_min_l, w_min_l) annotation (Line(points={{41,-4},
            {70,-4},{70,40},{110,40}},
                                   color={0,0,127}));
    connect(minWork4DOASCooling.w_min_sa, w_min_sa) annotation (Line(points={{41,-8},
            {110,-8}},                    color={0,0,127}));
    connect(minWork4DOASCooling.w_min_hw, w_min_hw) annotation (Line(points={{
            41,-12},{70,-12},{70,-40},{110,-40}}, color={0,0,127}));
    connect(minWork4DOASCooling.w_min_erw, w_min_erw) annotation (Line(points={
            {41,-16},{60,-16},{60,-74},{110,-74}}, color={0,0,127}));
      annotation (Line(points={{18,-14},{18,-10},{18,-10}}, color={0,0,127}),
              experiment(Tolerance=1e-06, __Dymola_Algorithm="Cvode"));
  end MinWork;

  model BaselineChiller "Baseline Chiller System"

    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.VolumeFlowRate V_w = 1e-2
        "volumetric flow rate of the chilled water";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of the heat exchanger";
      //parameter Real unit = 1.0 "number of parallel terminals";
      parameter Unit.Power W_fan = 10.0 "fan power";

      // Control Parameters
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Real t_lcmin = 0.1 "minimum leaving chilled water temperature in degC";
      parameter Unit.Efficiency deltaCOP = 1e-2;
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of the BCS in kW"
      annotation (Placement(transformation(extent={{100,60},{120,80}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of the BCS"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{60,20},{80,40}})));
    Chiller chiller(P_1 = P, P_4 = P, V_a = V_a, V_w = V_w, epsilon = epsilon,
      deltaQ = deltaQ, t_lcmin = t_lcmin)
      annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
    ChillerCOPFittingCurve chillerCOPFittingCurve(deltaCOP = deltaCOP, deltaQ = deltaQ)
      annotation (Placement(transformation(extent={{20,-40},{40,-20}})));


  equation

    connect(q_lat, dOASCoolingLoad.q_lat)
      annotation (Line(points={{-120,80},{-82,80},{-82,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-88,40},
            {-88,34.6},{-82,34.6}}, color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-100,0},{
            -100,30},{-82,30}}, color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3)
      annotation (Line(points={{-120,-80},{-82,-80},{-82,21}}, color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-92,-40},
            {-92,25.4},{-82,25.4}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-59,36},{58,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.w_set, chiller.w_set) annotation (Line(points={{-59,30},
            {-50,30},{-50,-38},{-42,-38}}, color={0,0,127}));
    connect(t_o, chiller.t_1) annotation (Line(points={{-120,40},{-88,40},{-88,-22},
            {-42,-22}}, color={0,0,127}));
    connect(t_i, chiller.t_4) annotation (Line(points={{-120,0},{-100,0},{-100,-27.4},
            {-42,-27.4}}, color={0,0,127}));
    connect(rh_o, chiller.rh_1) annotation (Line(points={{-120,-40},{-60,-40},{-60,
            -32.6},{-42,-32.6}}, color={0,0,127}));
    connect(chiller.q_c, chillerCOPFittingCurve.q_c) annotation (Line(points={{-19,-22},
            {-1.5,-22},{-1.5,-23},{18,-23}},          color={0,0,127}));
    connect(chiller.t_5, chillerCOPFittingCurve.t_l) annotation (Line(points={{-19,-30},
            {-0.5,-30},{-0.5,-30},{18,-30}},          color={0,0,127}));
    connect(chiller.t_wc, chillerCOPFittingCurve.t_e) annotation (Line(points={{-19,-38},
            {-0.5,-38},{-0.5,-37},{18,-37}},          color={0,0,127}));
    connect(chillerCOPFittingCurve.w_c, cOP.w_c) annotation (Line(points={{41,-34},
            {48,-34},{48,24},{58,24}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-59,36},{22,36},{22,
            70},{110,70}}, color={0,0,127}));
    connect(cOP.COP, COP)
      annotation (Line(points={{81,30},{110,30}}, color={0,0,127}));
  annotation (experiment(Tolerance=1e-06, __Dymola_Algorithm="Cvode"));
  end BaselineChiller;

  model ERWChiller "ERW + Chiller HVAC System"

    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.VolumeFlowRate V_w = 1e-2
        "volumetric flow rate of the chilled water";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of the heat exchanger";
      parameter Unit.Efficiency epsilon_s(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Efficiency epsilon_l(max = 1.0) = 0.75
        "effectiveness of latent heat recovery";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in ERW to supply air";
      //parameter Real unit = 1.0 "number of parallel terminals";
      parameter Unit.Power W_fan = 10.0 "fan power";
      parameter Unit.Pressure P_fan = 100
        "pressure loss across the HRW system";


      // Control Parameters
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Real t_lcmin = 0.1 "minimum leaving chilled water temperature in degC";
      parameter Unit.Efficiency deltaCOP = 1e-2 "number of parallel terminals";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,60},{120,80}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of ERW+Membrane(w/ condenser) system"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{64,20},{84,40}})));
    Chiller chiller(P_1 = P, P_4 = P, V_a = V_a, V_w = V_w, epsilon = epsilon,
      deltaQ = deltaQ, t_lcmin = t_lcmin)
      annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
    ChillerCOPFittingCurve chillerCOPFittingCurve(deltaCOP = deltaCOP)
      annotation (Placement(transformation(extent={{32,-40},{52,-20}})));

    HVACThermo.EnthalpyRecoveryWheel enthalpyRecoveryWheel(P = P, epsilon_s = epsilon_s,
      epsilon_l = epsilon_l, r = r, deltaT = deltaT, deltaw = deltaw, deltaP = deltaP,
      P_fan = P_fan, V_a = V_a)
      annotation (Placement(transformation(extent={{-62,-40},{-42,-20}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{24,-6},{44,14}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat)
      annotation (Line(points={{-120,80},{-82,80},{-82,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-88,40},
            {-88,34.6},{-82,34.6}}, color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-96,0},{-96,
            30},{-82,30}},      color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3)
      annotation (Line(points={{-120,-80},{-102,-80},{-102,-78},{-82,-78},{-82,21}},
                                                               color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-92,-40},
            {-92,25.4},{-82,25.4}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-59,36},{62,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.w_set, chiller.w_set) annotation (Line(points={{-59,30},
            {-32,30},{-32,-38},{-12,-38}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-59,36},{22,36},{22,
            70},{110,70}}, color={0,0,127}));
    connect(cOP.COP, COP)
      annotation (Line(points={{85,30},{110,30}}, color={0,0,127}));
    connect(t_o, enthalpyRecoveryWheel.t_1) annotation (Line(points={{-120,40},
            {-88,40},{-88,-22},{-64,-22}},
                                      color={0,0,127}));
    connect(t_i, enthalpyRecoveryWheel.t_2) annotation (Line(points={{-120,0},{
            -96,0},{-96,-27.2},{-64,-27.2}},
                                         color={0,0,127}));
    connect(rh_o, enthalpyRecoveryWheel.rh_1) annotation (Line(points={{-120,
            -40},{-92,-40},{-92,-32.8},{-64,-32.8}},
                                                color={0,0,127}));
    connect(rh_i, enthalpyRecoveryWheel.rh_2) annotation (Line(points={{-120,
            -80},{-82,-80},{-82,-38},{-64,-38}},
                                            color={0,0,127}));
    connect(t_i, chiller.t_4) annotation (Line(points={{-120,0},{-24,0},{-24,
            -27.4},{-12,-27.4}},
                          color={0,0,127}));
    connect(chiller.t_5, chillerCOPFittingCurve.t_l) annotation (Line(points={{11,-30},
            {20.5,-30},{20.5,-30},{30,-30}},             color={0,0,127}));
    connect(chiller.t_wc, chillerCOPFittingCurve.t_e) annotation (Line(points={{11,-38},
            {18,-38},{18,-37},{30,-37}},              color={0,0,127}));
    connect(chiller.q_c, chillerCOPFittingCurve.q_c) annotation (Line(points={{11,-22},
            {30,-22},{30,-23}},             color={0,0,127}));
    connect(enthalpyRecoveryWheel.rh_3, chiller.rh_1) annotation (Line(points={{-41,
            -34.4},{-12,-34.4},{-12,-32.6}},  color={0,0,127}));
    connect(enthalpyRecoveryWheel.w_p, add.u1) annotation (Line(points={{-41,-21},
            {-36.5,-21},{-36.5,10},{22,10}}, color={0,0,127}));
    connect(chillerCOPFittingCurve.w_c, add.u2) annotation (Line(points={{53,-34},
            {62,-34},{62,-10},{10,-10},{10,-2},{22,-2}}, color={0,0,127}));
    connect(add.y, cOP.w_c)
      annotation (Line(points={{45,4},{52,4},{52,24},{62,24}}, color={0,0,127}));
    connect(enthalpyRecoveryWheel.t_3, chiller.t_1) annotation (Line(points={{
            -41,-25.4},{-26.5,-25.4},{-26.5,-22},{-12,-22}}, color={0,0,127}));
  annotation (experiment(Tolerance=1e-06, __Dymola_Algorithm="Cvode"));
  end ERWChiller;



  model HRWChiller "HRW + Chiller HVAC System"

    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.VolumeFlowRate V_w = 1e-2
        "volumetric flow rate of the chilled water";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of the heat exchanger";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in ERW to supply air";
      //parameter Real unit = 1.0 "number of parallel terminals";
      parameter Unit.Power W_fan = 10.0 "fan power";
      parameter Unit.Pressure P_fan = 100
        "pressure loss across the HRW system";

      // Control Parameters
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Real t_lcmin = 0.1 "minimum leaving chilled water temperature in degC";
      parameter Unit.Efficiency deltaCOP = 1e-2 "number of parallel terminals";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,60},{120,80}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of ERW+Membrane(w/ condenser) system"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{64,20},{84,40}})));
    Chiller chiller(P_1 = P, P_4 = P, V_a = V_a, V_w = V_w, epsilon = epsilon,
      deltaQ = deltaQ, t_lcmin = t_lcmin)
      annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
    ChillerCOPFittingCurve chillerCOPFittingCurve(deltaCOP = deltaCOP)
      annotation (Placement(transformation(extent={{32,-40},{52,-20}})));

    HVACThermo.HeatRecoveryWheel heatRecoveryWheel(P = P, r = r, epsilon = epsilon,
      deltaT = deltaT, P_fan = P_fan, V_a = V_a)
      annotation (Placement(transformation(extent={{-62,-40},{-42,-20}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{20,-4},{40,16}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat)
      annotation (Line(points={{-120,80},{-82,80},{-82,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-88,40},
            {-88,34.6},{-82,34.6}}, color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-96,0},{-96,
            30},{-82,30}},      color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3)
      annotation (Line(points={{-120,-80},{-82,-80},{-82,21}}, color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-92,-40},
            {-92,25.4},{-82,25.4}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-59,36},{62,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.w_set, chiller.w_set) annotation (Line(points={{-59,30},
            {-32,30},{-32,-38},{-12,-38}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-59,36},{22,36},{22,
            70},{110,70}}, color={0,0,127}));
    connect(cOP.COP, COP)
      annotation (Line(points={{85,30},{110,30}}, color={0,0,127}));
    connect(t_i, chiller.t_4) annotation (Line(points={{-120,0},{-24,0},{-24,
            -27.4},{-12,-27.4}},
                          color={0,0,127}));
    connect(chiller.t_5, chillerCOPFittingCurve.t_l) annotation (Line(points={{11,-30},
            {20.5,-30},{20.5,-30},{30,-30}},             color={0,0,127}));
    connect(chiller.t_wc, chillerCOPFittingCurve.t_e) annotation (Line(points={{11,-38},
            {18,-38},{18,-37},{30,-37}},              color={0,0,127}));
    connect(chiller.q_c, chillerCOPFittingCurve.q_c) annotation (Line(points={{11,-22},
            {30,-22},{30,-23}},             color={0,0,127}));
    connect(t_o, heatRecoveryWheel.t_1) annotation (Line(points={{-120,40},{-88,40},
            {-88,-22},{-64,-22}}, color={0,0,127}));
    connect(t_i, heatRecoveryWheel.t_2) annotation (Line(points={{-120,0},{-96,0},
            {-96,-27.2},{-64,-27.2}}, color={0,0,127}));
    connect(rh_o, heatRecoveryWheel.rh_1) annotation (Line(points={{-120,-40},{-92,
            -40},{-92,-32.8},{-64,-32.8}}, color={0,0,127}));
    connect(rh_i, heatRecoveryWheel.rh_2) annotation (Line(points={{-120,-80},{-86,
            -80},{-86,-38},{-64,-38}}, color={0,0,127}));
    connect(heatRecoveryWheel.t_3, chiller.t_1) annotation (Line(points={{-41,-25.6},
            {-26.5,-25.6},{-26.5,-22},{-12,-22}},
                                                color={0,0,127}));
    connect(heatRecoveryWheel.rh_3, chiller.rh_1) annotation (Line(points={{-41,-34.4},
            {-26.5,-34.4},{-26.5,-32.6},{-12,-32.6}},
                                                    color={0,0,127}));
    connect(heatRecoveryWheel.w_p, add.u1) annotation (Line(points={{-41,-21.2},{-36.5,
            -21.2},{-36.5,12},{18,12}}, color={0,0,127}));
    connect(chillerCOPFittingCurve.w_c, add.u2) annotation (Line(points={{53,-34},
            {60,-34},{60,-10},{0,-10},{0,0},{18,0}}, color={0,0,127}));
    connect(add.y, cOP.w_c)
      annotation (Line(points={{41,6},{52,6},{52,24},{62,24}}, color={0,0,127}));
  annotation (experiment(Tolerance=1e-06, __Dymola_Algorithm="Cvode"));
  end HRWChiller;

  model DesiccantCooling
    "Desiccant Wheel + Heat Exchanger + GCHX HVAC System"
    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Efficiency epsilon_b(max = 1.0) = 0.80
        "effectiveness of ground heat exchanger";
      parameter Real theta_o = 1.1
        "dimensionless outlet temperature of DW";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air to supply air";
      parameter Unit.Efficiency COP_h(min = 0) = 4.8
        "COP of the heater";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.Power W_fan = 10.0 "fan power";
      parameter Unit.Diameter d_p = 2.5e-2 "diameter of pipe";
      parameter Unit.Height H = 30 "Depth of the borehole";
      parameter Real n = 1.0 "number of parallel boreholes";
      parameter Unit.Temperature t_g(min = 0) = 10.0
        "annual average ground temperature";
     parameter Real gCHXisOn = 1.0
        "whether groud heat exchanger is used in the analysis, 1 represent on";
      parameter Real r_R = 0.1
        "ratio of thermal resistance of borehole to the ground";
      parameter Unit.Pressure P_fan = 100.0 "pressure drop across the HRW";
      parameter Unit.Pressure P_dw = 120.0 "pressure drop across the DW";

      // Control Parameters
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the HX system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the desiccant wheel";
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of ERW+Membrane(w/ condenser) system"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{66,20},{86,40}})));
    HVACThermo.DesiccantWheel desiccantWheel(P_a = P, theta_o = theta_o, P_dw = P_dw,
      V_a = V_a, COP_h = COP_h, r = r, deltaw = deltaw, deltaWc = deltaWc)
      annotation (Placement(transformation(extent={{-60,-40},{-40,-20}})));
    HVACThermo.HeatRecoveryWheel heatRecoveryWheel(P = P, epsilon = epsilon, r = r,
      deltaT = deltaT, P_fan = P_fan, V_a = V_a)
      annotation (Placement(transformation(extent={{-16,-40},{4,-20}})));
    Modelica.Blocks.Sources.Constant q(k=0) "additional cooling load in kW"
      annotation (Placement(transformation(extent={{-8,-86},{12,-66}})));
    HVACThermo.TWetBulb tWetBulb(P = P)
      annotation (Placement(transformation(extent={{-40,60},{-20,80}})));
    Modelica.Blocks.Math.Add3 add3_1
      annotation (Placement(transformation(extent={{32,-10},{52,10}})));
    HVACThermo.GCHX gCHX(
      P=P,
      epsilon=epsilon,
      d_p=d_p,
      t_g=t_g,
      H=H,
      epsilon_b=epsilon_b,
      deltaT=deltaT,
      n=n,
      V_a=V_a,
      r_R=r_R,
      gCHXisOn = gCHXisOn)
      annotation (Placement(transformation(extent={{32,-40},{52,-20}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat) annotation (Line(points={{-120,80},{-82,
            80},{-82,39},{-42,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-86,40},
            {-86,34.6},{-42,34.6}}, color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-82,-40},
            {-82,25.4},{-42,25.4}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-19,36},{64,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-19,36},{0,36},{0,
            60},{110,60}}, color={0,0,127}));
    connect(cOP.COP, COP) annotation (Line(points={{87,30},{110,30}},
          color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-86,0},{-86,
            30},{-42,30}}, color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3) annotation (Line(points={{-120,-80},{
            -76,-80},{-76,21},{-42,21}},
                                color={0,0,127}));
    connect(dOASCoolingLoad.w_set, desiccantWheel.w_set) annotation (Line(points={
            {-19,30},{-12,30},{-12,-10},{-68,-10},{-68,-22},{-62,-22}}, color={0,0,
            127}));
    connect(t_o, desiccantWheel.t_1) annotation (Line(points={{-120,40},{-92,40},{
            -92,-26},{-62,-26}}, color={0,0,127}));
    connect(rh_o, desiccantWheel.rh_1) annotation (Line(points={{-120,-40},{-82,-40},
            {-82,-30},{-62,-30}}, color={0,0,127}));
    connect(desiccantWheel.t_2, heatRecoveryWheel.t_1) annotation (Line(points={{-39,
            -30},{-32,-30},{-32,-22},{-18,-22}}, color={0,0,127}));
    connect(desiccantWheel.rh_2, heatRecoveryWheel.rh_1) annotation (Line(points={
            {-39,-37},{-26.5,-37},{-26.5,-32.8},{-18,-32.8}}, color={0,0,127}));
    connect(rh_i, heatRecoveryWheel.rh_2) annotation (Line(points={{-120,-80},{-24,
            -80},{-24,-38},{-18,-38}}, color={0,0,127}));
    connect(t_i, heatRecoveryWheel.t_2) annotation (Line(points={{-120,0},{-28,0},
            {-28,-27.2},{-18,-27.2}}, color={0,0,127}));
    connect(heatRecoveryWheel.t_4, desiccantWheel.t_4) annotation (Line(points={{5,-30},
            {14,-30},{14,-56},{-72,-56},{-72,-34},{-62,-34}},          color={0,
            0,127}));
    connect(heatRecoveryWheel.rh_4, desiccantWheel.rh_4) annotation (Line(
          points={{5,-39},{10,-39},{10,-50},{-68,-50},{-68,-38},{-62,-38}},
          color={0,0,127}));
    connect(t_o, tWetBulb.t) annotation (Line(points={{-120,40},{-92,40},{-92,75},
            {-42,75}}, color={0,0,127}));
    connect(rh_o, tWetBulb.rh) annotation (Line(points={{-120,-40},{-78,-40},{-78,
            65},{-42,65}}, color={0,0,127}));
    connect(desiccantWheel.w_c, add3_1.u1) annotation (Line(points={{-39,-23},{-34,
            -23},{-34,8},{30,8}}, color={0,0,127}));
    connect(heatRecoveryWheel.w_p, add3_1.u2) annotation (Line(points={{5,-21.2},{
            8,-21.2},{8,0},{30,0}}, color={0,0,127}));
    connect(add3_1.y, cOP.w_c)
      annotation (Line(points={{53,0},{58,0},{58,24},{64,24}}, color={0,0,127}));
    connect(q.y, gCHX.q_c) annotation (Line(points={{13,-76},{20,-76},{20,-22},{30,
            -22}}, color={0,0,127}));
    connect(heatRecoveryWheel.t_3, gCHX.t_1) annotation (Line(points={{5,-25.6},{17.5,
            -25.6},{17.5,-26},{30,-26}}, color={0,0,127}));
    connect(t_i, gCHX.t_2) annotation (Line(points={{-120,0},{0,0},{0,-10},{24,-10},
            {24,-30},{30,-30}}, color={0,0,127}));
    connect(heatRecoveryWheel.rh_3, gCHX.rh_1) annotation (Line(points={{5,-34.4},
            {17.5,-34.4},{17.5,-34},{30,-34}}, color={0,0,127}));
    connect(tWetBulb.t_wb, gCHX.t_wb) annotation (Line(points={{-19,70},{16,70},{16,
            -38},{30,-38}}, color={0,0,127}));
    connect(gCHX.w_c, add3_1.u3) annotation (Line(points={{53,-30},{60,-30},{60,-14},
            {30,-14},{30,-8}}, color={0,0,127}));
    annotation (Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Desiccant Wheel + Heat Exchanger + Ground-Coupled Heat Exchanger with Condenser Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of a coupled system with an upstream desiccant wheel and a downstream heat exchanger to reuse the cool exhaust air, and a downstream ground-coupled heat exchanger (or an alternative chiller when ground temperature is high) to provide additional sensible cooling.</span></p>
<p><br><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end DesiccantCooling;

  model MembraneCooling
    "ERW + Membrane (w/ condenser) + GCHX HVAC System"
    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.Efficiency epsilon_s(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Efficiency epsilon_l(max = 1.0) = 0.75
        "effectiveness of latent heat recovery";
      parameter Unit.Efficiency eta_is(max = 1.0) = 0.9
        "isentropic efficiency of the vacuum pump";
      parameter Unit.Pressure P_stage = 2000 "pressure at the exit of stage pump";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in ERW to supply air";
      parameter Unit.Power W_fan = 10.0 "fan power";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Diameter d_p = 2.5e-2 "diameter of pipe";
      parameter Unit.Height H = 30 "depth of the borehole";
      parameter Unit.Efficiency epsilon_b(max = 1.0) = 0.80
        "effectiveness of ground heat exchanger";
      parameter Real n = 1.0 "number of parallel boreholes";
      parameter Unit.Temperature t_g(min = 0) = 10.0
        "annual average ground temperature";
      parameter Real gCHXisOn = 1.0
        "whether groud heat exchanger is used in the analysis, 1 represent on";
      parameter Real r_R = 0.1
        "ratio of thermal resistance of borehole to the ground";
      parameter Unit.Area A_m = 1.0
        "Membrane area per channel in m2";
      parameter Real n_m = 100.0 "number of membrane stacks";
      parameter Unit.Length h = 0.01 "Channel height in m";
      parameter Real p_w = 6.2e-7 "water permeance of the membrane in mol/m2 s Pa";
      parameter Unit.Pressure P_fan = 100
        "pressure loss across the ERV system";


      // Control Parameters
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";

    MembraneWithCondenser membraneWithCondenser(P_a = P, eta_is = eta_is,
      V_a = V_a, P_stage = P_stage, deltaw_set = deltaw, deltaP_stage = deltaP,
      deltaWc = deltaWc, A_m = A_m, n_m = n_m, h = h, p_w = p_w)
      annotation (Placement(transformation(extent={{-6,-40},{14,-20}})));
    EnthalpyRecoveryWheel enthalpyRecoveryWheel(P = P, epsilon_s = epsilon_s,
      epsilon_l = epsilon_l, r = r, deltaT = deltaT, deltaw = deltaw,
      deltaP = deltaP, V_a = V_a, P_fan = P_fan)
      annotation (Placement(transformation(extent={{-46,-40},{-26,-20}})));
    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of ERW+Membrane(w/ condenser) system"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{70,20},{90,40}})));
    HVACThermo.GCHX gCHX(P = P, epsilon = epsilon, d_p = d_p, t_g = t_g, H = H,
      epsilon_b = epsilon_b, deltaT = deltaT, n = n, V_a = V_a, r_R = r_R,
      gCHXisOn = gCHXisOn)
      annotation (Placement(transformation(extent={{32,-40},{52,-20}})));
    HVACThermo.TWetBulb tWetBulb(P = P)
      annotation (Placement(transformation(extent={{-10,-76},{10,-56}})));
    Modelica.Blocks.Math.Add3 add3_1
      annotation (Placement(transformation(extent={{24,8},{42,26}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat) annotation (Line(points={{-120,80},{-82,
            80},{-82,39},{-42,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-86,40},
            {-86,34.6},{-42,34.6}}, color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-82,-40},
            {-82,25.4},{-42,25.4}}, color={0,0,127}));
    connect(t_o, enthalpyRecoveryWheel.t_1) annotation (Line(points={{-120,40},{-54,
            40},{-54,-22},{-48,-22}}, color={0,0,127}));
    connect(t_i, enthalpyRecoveryWheel.t_2) annotation (Line(points={{-120,0},{-60,
            0},{-60,-27.2},{-48,-27.2}}, color={0,0,127}));
    connect(rh_o, enthalpyRecoveryWheel.rh_1) annotation (Line(points={{-120,-40},
            {-60,-40},{-60,-32.8},{-48,-32.8}}, color={0,0,127}));
    connect(rh_i, enthalpyRecoveryWheel.rh_2) annotation (Line(points={{-120,-80},
            {-54,-80},{-54,-38},{-48,-38}}, color={0,0,127}));
    connect(dOASCoolingLoad.w_set, membraneWithCondenser.w_set) annotation (Line(
          points={{-19,30},{-16,30},{-16,-22},{-8,-22}},
                                                     color={0,0,127}));
    connect(enthalpyRecoveryWheel.t_3, membraneWithCondenser.t_1) annotation (
        Line(points={{-25,-25.4},{-20.5,-25.4},{-20.5,-30},{-8,-30}},
                                                                color={0,0,127}));
    connect(enthalpyRecoveryWheel.rh_3, membraneWithCondenser.rh_1) annotation (
        Line(points={{-25,-34.4},{-18,-34.4},{-18,-38},{-8,-38}},
                                                                color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-19,36},{68,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-19,36},{0,36},{0,
            60},{110,60}}, color={0,0,127}));
    connect(cOP.COP, COP) annotation (Line(points={{91,30},{110,30}},
          color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-86,0},{-86,
            30},{-42,30}}, color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3) annotation (Line(points={{-120,-80},{-72,-80},
            {-72,21},{-42,21}}, color={0,0,127}));
    connect(membraneWithCondenser.rh_s, gCHX.rh_1) annotation (Line(points={{15,-38},
            {24.5,-38},{24.5,-34},{30,-34}},     color={0,0,127}));
    connect(membraneWithCondenser.t_s, gCHX.t_1) annotation (Line(points={{15,-33},
            {22.5,-33},{22.5,-26},{30,-26}},   color={0,0,127}));
    connect(t_i, gCHX.t_2) annotation (Line(points={{-120,0},{26,0},{26,-30},{30,-30}},
                     color={0,0,127}));
    connect(membraneWithCondenser.q_cond, gCHX.q_c)
      annotation (Line(points={{15,-22},{30,-22}}, color={0,0,127}));
    connect(t_o, tWetBulb.t) annotation (Line(points={{-120,40},{-66,40},{-66,-61},
            {-12,-61}}, color={0,0,127}));
    connect(rh_o, tWetBulb.rh) annotation (Line(points={{-120,-40},{-82,-40},{-82,
            -71},{-12,-71}}, color={0,0,127}));
    connect(tWetBulb.t_wb, gCHX.t_wb) annotation (Line(points={{11,-66},{24,-66},{
            24,-42},{30,-42},{30,-38}}, color={0,0,127}));
    connect(add3_1.y, cOP.w_c) annotation (Line(points={{42.9,17},{55.45,17},{55.45,
            24},{68,24}}, color={0,0,127}));
    connect(enthalpyRecoveryWheel.w_p, add3_1.u1) annotation (Line(points={{-25,-21},
            {-20,-21},{-20,12},{0,12},{0,24.2},{22.2,24.2}}, color={0,0,127}));
    connect(membraneWithCondenser.w_c, add3_1.u2) annotation (Line(points={{15,-27},
            {20,-27},{20,-8},{6,-8},{6,17},{22.2,17}}, color={0,0,127}));
    connect(gCHX.w_c, add3_1.u3) annotation (Line(points={{53,-30},{60,-30},{60,4},
            {16,4},{16,9.8},{22.2,9.8}}, color={0,0,127}));
    annotation (Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Enthalpy Recovery Wheel + Membrane Dehumidifier + Ground-Coupled Heat Exchanger with Condenser Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of a coupled system with an upstream ERW and a downstream membrane dehumidifier with a condenser recycling excess water vapor in the air, and a downstream ground-coupled heat exchanger.</span></p>
<p><br><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end MembraneCooling;

  model MembraneCoolingWithMemb
     "Membrane (w/ Membrane) + HRW + GCHX HVAC System"
    import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.Efficiency eta_is(max = 1.0) = 0.9
        "isentropic efficiency of the vacuum pump";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in ERW to supply air";
      parameter Unit.Power W_fan = 10.0 "fan power";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.Diameter d_p = 2.5e-2 "diameter of pipe";
      parameter Unit.Height H = 30 "depth of the borehole";
      parameter Unit.Efficiency epsilon_b(max = 1.0) = 0.80
        "effectiveness of ground heat exchanger";
      parameter Real n = 1.0 "number of parallel boreholes";
      parameter Unit.Temperature t_g(min = 0) = 10.0
        "annual average ground temperature";
      parameter Real gCHXisOn = 1.0
        "whether groud heat exchanger is used in the analysis, 1 represent on";
      parameter Real r_R = 0.1
        "ratio of thermal resistance of borehole to the ground";
      parameter Unit.Area A_m = 1.0
        "Membrane area per channel in m2";
      parameter Real n_m = 100.0 "number of membrane stacks";
      parameter Unit.Length h = 0.01 "Channel height in m";
      parameter Real p_w = 6.2e-7 "water permeance of the membrane in mol/m2 s Pa";
      parameter Unit.Pressure P_fan = 100
        "pressure loss across the HRW";
      parameter Unit.Efficiency COP_h = 4.8 "COP of heat pump";


      // Control Parameters
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the ERV system";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the ERV system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of ERW+Membrane(w/ condenser) system"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{70,20},{90,40}})));
    HVACThermo.GCHX gCHX(P = P, epsilon = epsilon, d_p = d_p, t_g = t_g, H = H,
      epsilon_b = epsilon_b, deltaT = deltaT, n = n, V_a = V_a, r_R = r_R,
      gCHXisOn = gCHXisOn)
      annotation (Placement(transformation(extent={{32,-40},{52,-20}})));
    HVACThermo.TWetBulb tWetBulb(P = P)
      annotation (Placement(transformation(extent={{-10,-76},{10,-56}})));
    Modelica.Blocks.Math.Add3 add3_1
      annotation (Placement(transformation(extent={{24,8},{42,26}})));
    HVACThermo.MembraneWithMembrane membraneWithMembrane(P_a = P, V_a = V_a,
      deltaw_set = deltaw, deltaWc = deltaWc, A_m = A_m, n_m = n_m, h = h,
      p_w = p_w, COP_h = COP_h, eta_is = eta_is, deltaP_stage = deltaP)
      annotation (Placement(transformation(extent={{-48,-40},{-28,-20}})));
    HVACThermo.HeatRecoveryWheel heatRecoveryWheel(P = P, r = r, epsilon = epsilon,
      deltaT = deltaT, P_fan = P_fan, V_a = V_a)
      annotation (Placement(transformation(extent={{-8,-40},{12,-20}})));
    Modelica.Blocks.Sources.Constant q(k=0) "additional cooling load in kW"
      annotation (Placement(transformation(extent={{40,-76},{60,-56}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat) annotation (Line(points={{-120,80},{-82,
            80},{-82,39},{-42,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-86,40},
            {-86,34.6},{-42,34.6}}, color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-82,-40},
            {-82,25.4},{-42,25.4}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-19,36},{68,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-19,36},{0,36},{0,
            60},{110,60}}, color={0,0,127}));
    connect(cOP.COP, COP) annotation (Line(points={{91,30},{110,30}},
          color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-86,0},{-86,
            30},{-42,30}}, color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3) annotation (Line(points={{-120,-80},{-72,-80},
            {-72,21},{-42,21}}, color={0,0,127}));
    connect(t_i, gCHX.t_2) annotation (Line(points={{-120,0},{26,0},{26,-30},{30,-30}},
                     color={0,0,127}));
    connect(t_o, tWetBulb.t) annotation (Line(points={{-120,40},{-66,40},{-66,-61},
            {-12,-61}}, color={0,0,127}));
    connect(rh_o, tWetBulb.rh) annotation (Line(points={{-120,-40},{-82,-40},{-82,
            -71},{-12,-71}}, color={0,0,127}));
    connect(add3_1.y, cOP.w_c) annotation (Line(points={{42.9,17},{55.45,17},{55.45,
            24},{68,24}}, color={0,0,127}));
    connect(gCHX.w_c, add3_1.u3) annotation (Line(points={{53,-30},{60,-30},{60,4},
            {16,4},{16,9.8},{22.2,9.8}}, color={0,0,127}));
    connect(dOASCoolingLoad.w_set, membraneWithMembrane.w_set) annotation (Line(
          points={{-19,30},{-14,30},{-14,-12},{-60,-12},{-60,-22},{-50,-22}},
          color={0,0,127}));
    connect(tWetBulb.t_wb, gCHX.t_wb) annotation (Line(points={{11,-66},{26,-66},{
            26,-38},{30,-38}}, color={0,0,127}));
    connect(t_o, membraneWithMembrane.t_1) annotation (Line(points={{-120,40},{-66,
            40},{-66,-26},{-50,-26}}, color={0,0,127}));
    connect(rh_o, membraneWithMembrane.rh_1) annotation (Line(points={{-120,-40},{
            -82,-40},{-82,-30},{-50,-30}}, color={0,0,127}));
    connect(membraneWithMembrane.w_c, add3_1.u1) annotation (Line(points={{-27,-22},
            {-20,-22},{-20,14},{0,14},{0,24.2},{22.2,24.2}}, color={0,0,127}));
    connect(heatRecoveryWheel.w_p, add3_1.u2) annotation (Line(points={{13,-21.2},
            {20,-21.2},{20,-10},{8,-10},{8,17},{22.2,17}}, color={0,0,127}));
    connect(heatRecoveryWheel.t_3, gCHX.t_1) annotation (Line(points={{13,-25.6},{
            21.5,-25.6},{21.5,-26},{30,-26}}, color={0,0,127}));
    connect(heatRecoveryWheel.rh_3, gCHX.rh_1) annotation (Line(points={{13,-34.4},
            {21.5,-34.4},{21.5,-34},{30,-34}}, color={0,0,127}));
    connect(heatRecoveryWheel.t_4, membraneWithMembrane.t_5) annotation (Line(
          points={{13,-30},{20,-30},{20,-50},{-60,-50},{-60,-34},{-50,-34}},
          color={0,0,127}));
    connect(heatRecoveryWheel.rh_4, membraneWithMembrane.rh_5) annotation (Line(
          points={{13,-39},{18,-39},{18,-46},{-56,-46},{-56,-38},{-50,-38}},
          color={0,0,127}));
    connect(membraneWithMembrane.t_s, heatRecoveryWheel.t_1) annotation (Line(
          points={{-27,-30},{-16,-30},{-16,-22},{-10,-22}}, color={0,0,127}));
    connect(t_i, heatRecoveryWheel.t_2) annotation (Line(points={{-120,0},{-24,0},
            {-24,-27.2},{-10,-27.2}}, color={0,0,127}));
    connect(membraneWithMembrane.rh_s, heatRecoveryWheel.rh_1) annotation (Line(
          points={{-27,-38},{-20,-38},{-20,-32.8},{-10,-32.8}}, color={0,0,127}));
    connect(rh_i, heatRecoveryWheel.rh_2) annotation (Line(points={{-120,-80},{-72,
            -80},{-72,-56},{-16,-56},{-16,-38},{-10,-38}}, color={0,0,127}));
    connect(q.y, gCHX.q_c) annotation (Line(points={{61,-66},{80,-66},{80,-12},{24,
            -12},{24,-22},{30,-22}}, color={0,0,127}));
    annotation (Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Enthalpy Recovery Wheel + Membrane Dehumidifier + Ground-Coupled Heat Exchanger with Condenser Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of a coupled system with an upstream ERW and a downstream membrane dehumidifier with a condenser recycling excess water vapor in the air, and a downstream ground-coupled heat exchanger.</span></p>
<p><br><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end MembraneCoolingWithMemb;

  model MembraneCooling_SG
     "Membrane (w/ sweeping gas) + Chiller System"
      import HVACThermo;

      // Import package
      package Unit = Modelica.SIunits;

      // Parameter of the model
      parameter Unit.Pressure P = 101325 "outdoor air pressure";
      parameter Unit.Efficiency epsilon(max = 1.0) = 0.75
        "effectiveness of sensible heat recovery";
      parameter Unit.VolumeFlowRate V_a = 1.0
        "volumetric flow rate of the moist air";
      parameter Unit.VolumeFlowRate V_w = 1e-2
        "volumetric flow rate of the chilled water";
      parameter Real r(min = 0.0, max = 1.0) = 1.0
        "ratio of mass flow rate of exhaust air in ERW to supply air";
      parameter Unit.Power W_fan = 10.0 "fan power";
      parameter Unit.Area A_m = 1.0
        "Membrane area per channel in m2";
      parameter Real n_m = 100.0 "number of membrane stacks";
      parameter Unit.Length h = 0.01 "Channel height in m";
      parameter Real p_w = 6.2e-7 "water permeance of the membrane in mol/m2 s Pa";
      parameter Unit.Diameter d_p = 2.5e-2 "diameter of pipe";
      parameter Unit.Height H = 30 "Depth of the borehole";
      parameter Real n = 1.0 "number of parallel boreholes";
      parameter Unit.Temperature t_g(min = 0) = 10.0
        "annual average ground temperature";
      parameter Real r_R = 0.1
        "ratio of thermal resistance of borehole to the ground";
      parameter Unit.Efficiency COP_h(min = 0) = 4.8
        "COP of the heat pump";
      parameter Unit.Efficiency epsilon_b(max = 1.0) = 0.80
        "effectiveness of ground heat exchanger";
      parameter Unit.Efficiency eta_is(max = 1.0) = 0.9
        "isentropic efficiency of the vacuum pump";

      // Control Parameters
      parameter Real t_lcmin = 0.1 "minimum leaving chilled water temperature in degC";
      parameter Real deltaw(min = 0) = 1e-5
        "humidity threshold control on/off for the membrane system";
      parameter Real deltaP(min = 0) = 1
        "pressure threshold control on/off for the ERV system";
      parameter Unit.Efficiency deltaCOP = 1e-2;
      parameter Real phi_min = 0.1 "humidity control on supply air";
      parameter Unit.Power deltaQ = 1.0 "threshold control on total cooling load";
      parameter Unit.Power deltaWc = 1.0
        "threshold control on minimum cooling work";
      parameter Unit.Temperature deltaT(min = 0) = 0.1
        "temperature threshold control on/off for the HX system";

    Modelica.Blocks.Interfaces.RealInput q_lat "latent heat removal"
      annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
    Modelica.Blocks.Interfaces.RealInput t_o "outdoor air temperature"
      annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
    Modelica.Blocks.Interfaces.RealInput t_i "indoor air temperature"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput rh_o
      "outdoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
    Modelica.Blocks.Interfaces.RealInput rh_i "indoor air relative humidity in %"
      annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
    Modelica.Blocks.Interfaces.RealOutput q_c
      "cooling load of ERW+Membrane(w/ condenser) system in kW"
      annotation (Placement(transformation(extent={{100,50},{120,70}})));
    Modelica.Blocks.Interfaces.RealOutput COP
      "COP of ERW+Membrane(w/ condenser) system"
      annotation (Placement(transformation(extent={{100,20},{120,40}})));
    DOASCoolingLoad dOASCoolingLoad(P_1 = P, P_3 = P, V_a = V_a, deltaQ = deltaQ,
      phi_min = phi_min)
      annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
    HVACThermo.COP cOP(deltaQ = deltaQ, deltaWc = deltaWc, W_fan = W_fan)
      annotation (Placement(transformation(extent={{70,20},{90,40}})));
    HVACThermo.MembraneWithSweepGas membraneWithSweepGas(P_a = P, V_a = V_a, r = r,
      deltaw_set = deltaw, deltaWc = deltaWc, A_m = A_m, n_m = n_m, h = h,
      p_w = p_w, COP_h = COP_h, eta_is = eta_is)
      annotation (Placement(transformation(extent={{-40,-40},{-20,-20}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{20,0},{40,20}})));
    HVACThermo.GCHX gCHX(P = P, epsilon = epsilon, d_p = d_p, t_g = t_g, H = H,
      epsilon_b = epsilon_b, deltaT = deltaT, n = n, V_a = V_a, r_R = r_R)
      annotation (Placement(transformation(extent={{4,-40},{24,-20}})));
    Modelica.Blocks.Sources.Constant q(k=0) "additional cooling load in kW"
      annotation (Placement(transformation(extent={{20,-80},{40,-60}})));
    HVACThermo.TWetBulb tWetBulb(P=P)
      annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
  equation

    connect(q_lat, dOASCoolingLoad.q_lat) annotation (Line(points={{-120,80},{-82,
            80},{-82,39},{-42,39}}, color={0,0,127}));
    connect(t_o, dOASCoolingLoad.t_1) annotation (Line(points={{-120,40},{-86,40},
            {-86,34.6},{-42,34.6}}, color={0,0,127}));
    connect(rh_o, dOASCoolingLoad.rh_1) annotation (Line(points={{-120,-40},{-82,-40},
            {-82,25.4},{-42,25.4}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, cOP.q_c)
      annotation (Line(points={{-19,36},{68,36}}, color={0,0,127}));
    connect(dOASCoolingLoad.q_c, q_c) annotation (Line(points={{-19,36},{0,36},{0,
            60},{110,60}}, color={0,0,127}));
    connect(cOP.COP, COP) annotation (Line(points={{91,30},{110,30}},
          color={0,0,127}));
    connect(t_i, dOASCoolingLoad.t_3) annotation (Line(points={{-120,0},{-86,0},{-86,
            30},{-42,30}}, color={0,0,127}));
    connect(rh_i, dOASCoolingLoad.rh_3) annotation (Line(points={{-120,-80},{-72,-80},
            {-72,21},{-42,21}}, color={0,0,127}));
    connect(dOASCoolingLoad.w_set, membraneWithSweepGas.w_set) annotation (Line(
          points={{-19,30},{-14,30},{-14,-6},{-60,-6},{-60,-22},{-42,-22}}, color=
           {0,0,127}));
    connect(t_o, membraneWithSweepGas.t_1) annotation (Line(points={{-120,40},{-92,
            40},{-92,-26},{-42,-26}}, color={0,0,127}));
    connect(rh_o, membraneWithSweepGas.rh_1) annotation (Line(points={{-120,-40},{
            -82,-40},{-82,-30},{-42,-30}}, color={0,0,127}));
    connect(t_i, membraneWithSweepGas.t_3) annotation (Line(points={{-120,0},{-86,
            0},{-86,-34},{-42,-34}}, color={0,0,127}));
    connect(rh_i, membraneWithSweepGas.rh_3) annotation (Line(points={{-120,-80},{
            -72,-80},{-72,-38},{-42,-38}}, color={0,0,127}));
    connect(membraneWithSweepGas.w_c, add.u1) annotation (Line(points={{-19,-22},{
            -14,-22},{-14,-10},{0,-10},{0,16},{18,16}}, color={0,0,127}));
    connect(add.y, cOP.w_c) annotation (Line(points={{41,10},{54,10},{54,24},{68,24}},
          color={0,0,127}));
    connect(gCHX.w_c, add.u2) annotation (Line(points={{25,-30},{30,-30},{30,-10},
            {8,-10},{8,4},{18,4}}, color={0,0,127}));
    connect(q.y, gCHX.q_c) annotation (Line(points={{41,-70},{48,-70},{48,-52},{-6,
            -52},{-6,-22},{2,-22}}, color={0,0,127}));
    connect(t_i, gCHX.t_2) annotation (Line(points={{-120,0},{-10,0},{-10,-30},{2,
            -30}}, color={0,0,127}));
    connect(membraneWithSweepGas.t_s, gCHX.t_1) annotation (Line(points={{-19,-30},
            {-12,-30},{-12,-26},{2,-26}}, color={0,0,127}));
    connect(gCHX.t_wb, tWetBulb.t_wb) annotation (Line(points={{2,-38},{-10,-38},{
            -10,-70},{-19,-70}}, color={0,0,127}));
    connect(t_o, tWetBulb.t) annotation (Line(points={{-120,40},{-92,40},{-92,-65},
            {-42,-65}}, color={0,0,127}));
    connect(rh_o, tWetBulb.rh) annotation (Line(points={{-120,-40},{-82,-40},{-82,
            -75},{-42,-75}}, color={0,0,127}));
    connect(membraneWithSweepGas.rh_s, gCHX.rh_1) annotation (Line(points={{-19,-38},
            {-14,-38},{-14,-34},{2,-34}}, color={0,0,127}));
    annotation (Documentation(info="<html>
<p><b><span style=\"font-family: Arial;\">Membrane Dehumidifier With Sweeping Gas + Chiller Model</span></b></p>
<p><span style=\"font-family: Arial;\">The model describes the thermodynamic performance of a coupled system, featuring a membrane dehumidifier with sweeping gas and a downstream chiller to dehumidify the exhuast air.</span></p>
<p><br><span style=\"font-family: Arial;\">The model is constructed by Tianyi Chen, PhD candidate at MIT under the supervision of Prof. Leslie Norford.</span></p>
</html>"));
  end MembraneCooling_SG;

  model ChillerCOP_sensitivity

    // Import the package
    package Unit = Modelica.SIunits;
    package Data = Buildings.Fluid.Chillers.Data.ElectricEIR;
    package Func = Buildings.Utilities.Math.Functions;

    // inputs
    Modelica.Blocks.Interfaces.RealInput plr(min=0)
      "part load ratio"
        annotation (Placement(transformation(extent={{-140,50},{-100,90}})));
    Modelica.Blocks.Interfaces.RealInput t_l(min=0)
      "leaving chilled water temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealInput t_e(min=0)
      "entering condenser water temperature in degree C"
        annotation (Placement(transformation(extent={{-140,-90},{-100,-50}})));

    // Output
    Modelica.Blocks.Interfaces.RealOutput COP = eta
      "COP of the chiller"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));

    // Chiller model type
    Data.ElectricEIRChiller_Trane_CVHF_2567kW_11_77COP_VSD per;

    // Conversion of inputs to physical parameters
    Unit.Temperature T_lmin(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TEvaLvgMin)
      "lower limit on leaving chilled water temperature";
    Unit.Temperature T_lmax(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TEvaLvgMax)
      "upper limit on leaving chilled water temperature";
    Unit.Temperature T_emin(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TConEntMin)
      "lower limit on entering condenser water temperature";
    Unit.Temperature T_emax(displayUnit = "degC")=
      Unit.Conversions.to_degC(per.TConEntMax)
      "upper limit on entering condenser water temperature";

    // Regression Model for the Chiller
    //Unit.Temperature Tl = min(T_lmax, max(T_lmin, t_l))
    //  "actual leaving chilled water temperature";
    //Unit.Temperature Te = min(T_emax, max(T_emin, t_e))
    //  "actual entering condenser water temperature";

    Unit.Temperature Tl = t_l
      "actual leaving chilled water temperature";
    Unit.Temperature Te = t_e
      "actual entering condenser water temperature";



    // Chiller chiller capacity fraction
    Unit.Efficiency capFunT;

    // Chiller energy input ratio
    Unit.Efficiency EIRFunT;

    // Chiller energy input ratio
    Unit.Efficiency EIRFunPLR;

    // Chiller Part Load Ratio (PLR)
    Unit.Efficiency PLR;

    // COP of the chiller
    Unit.Efficiency eta;

  equation
    // Chiller chiller capacity fraction biquadratic curve
    capFunT = Func.smoothMax(
      x1 = 1E-6,
      x2 = Func.biquadratic(a = per.capFunT, x1 = Tl, x2 = Te),
      deltaX = 1E-7);

    // Chiller energy input ratio biquadratic curve
    EIRFunT = Func.biquadratic(a = per.EIRFunT, x1 = Tl, x2 = Te);

    // Chiller energy input ratio quadratic curve
    PLR = min(per.PLRMax, max(per.PLRMin, plr));
    EIRFunPLR = per.EIRFunPLR[1] + per.EIRFunPLR[2] * PLR +
      per.EIRFunPLR[3] * PLR^2;

    // Calculate COP of the Chiller
    eta = per.COP_nominal / EIRFunT / EIRFunPLR * PLR;

    annotation (Icon(graphics={
        Text(
          extent={{-160,140},{168,104}},
          lineColor={28,108,200},
          textString="%name")}));
  end ChillerCOP_sensitivity;

  function exergy_moistAir
    "thermodynamic exergy as function of temperature and humidity ratio"

    extends Modelica.Icons.Function;
    package Unit = Modelica.SIunits;

    input Unit.Temperature T "thermodynamic state temperature in K";
    input Real w "thermodynamic state humidity ratio";
    input Unit.Temperature T_0 "environment temperature in K";
    input Unit.Temperature w_0 "environment humidity ratio";
    output Unit.SpecificEnergy ksi "Specific enthalpy";

    parameter Unit.SpecificHeatCapacity Cp_a = 1003 "specific heat capacity of moist air";
    parameter Unit.SpecificHeatCapacity Cp_w = 1872 "specific heat capacity of water vapor";
    parameter Unit.SpecificHeatCapacity R_a = 287;

  algorithm
    ksi := (Cp_a + w*Cp_w) * T_0 * (T/T_0 - 1 - Modelica.Math.log(T/T_0))
            + R_a * T_0 * ((1 + 1.608*w) * Modelica.Math.log((1 + 1.608*w_0)/(1 + 1.608*w))
                          + 1.608*w * Modelica.Math.log(w/w_0));

    annotation (Inline=true);
  end exergy_moistAir;
  annotation (uses(Modelica(version="3.2.2"), Buildings(version="5.0.1")));
end HVACThermo;
