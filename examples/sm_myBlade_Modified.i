Begin Sierra Job

  define direction x with vector 1.0 0.0 0.0
  define direction y with vector 0.0 1.0 0.0
  define direction z with vector 0.0 0.0 1.0

  define point ptO with coordinates 0.0 0.0 0.0
  define point ptZ with coordinates 0.0 0.0 1.0
  define point ptX with coordinates 1.0 0.0 0.0

  define coordinate system sysR rectangular with point ptO point ptZ point ptX

   ############### Function Definitions ########################
   # { tf = 1.0 } 
   # { umax = 1.0e0 } 
   begin function u3
      type = analytic
      expression variable: t = GLOBAL time
      evaluate expression is "(t*{umax}/{tf})"
   end function u3

   begin function ramp
      type = analytic
      expression variable: t = GLOBAL time
      evaluate expression is "(t*{1}/{tf})"
   end function ramp

   begin function apply_force_x
     type = analytic
      expression variable: t = GLOBAL time
      expression variable: f = NODAL node_forces_var(x)
      evaluate expression is "(t*{1}/{tf}*f)" 
   end

     begin function apply_force_y
     type = analytic
      expression variable: t = GLOBAL time
      expression variable: f = NODAL node_forces_var(y)
      evaluate expression is "(t*{1}/{tf}*f)" 
   end

   begin function apply_force_z
     type = analytic
      expression variable: t = GLOBAL time
      expression variable: f = NODAL node_forces_var(z)
      evaluate expression is "(t*{1}/{tf}*f)" 
   end 

begin property specification for material Adhesive
   DENSITY      = 1100.0
    begin parameters for model elastic_orthotropic
        youngs modulus = 70e9
        poissons ratio  = 0.33
        E11          = 4560000.0
        E22          = 4560000.0
        E33          = 4560000.0
        NU12         = 0.3
        NU13         = 0.3
        NU23         = 0.3
        G12          = 1450000.0
        G13          = 1450000.0
        G23          = 1450000.0

        # Coordinate system
        COORDINATE SYSTEM = sysR
    end parameters for model elastic_orthotropic
end property specification for material Adhesive


begin property specification for material glass_triax
   DENSITY      = 1940.0
    begin parameters for model elastic_orthotropic
        youngs modulus = 70e9
        poissons ratio  = 0.33
        E11          = 28211400000.0
        E22          = 16238800000.0
        E33          = 15835500000.0
        NU12         = 0.497511
        NU13         = 0.18091
        NU23         = 0.27481
        G12          = 8248220000.0
        G13          = 3491240000.0
        G23          = 3491240000.0

        # Coordinate system
        COORDINATE SYSTEM = sysR
    end parameters for model elastic_orthotropic
end property specification for material glass_triax


begin property specification for material medium_density_foam
   DENSITY      = 130.0
    begin parameters for model elastic_orthotropic
        youngs modulus = 70e9
        poissons ratio  = 0.33
        E11          = 142500000.0
        E22          = 142500000.0
        E33          = 142500000.0
        NU12         = 0.3194
        NU13         = 0.3194
        NU23         = 0.3194
        G12          = 54001819.0
        G13          = 54001819.0
        G23          = 54001819.0

        # Coordinate system
        COORDINATE SYSTEM = sysR
    end parameters for model elastic_orthotropic
end property specification for material medium_density_foam


begin property specification for material glass_uni
   DENSITY      = 1940.0
    begin parameters for model elastic_orthotropic
        youngs modulus = 70e9
        poissons ratio  = 0.33
        E11          = 43700000000.0
        E22          = 16500000000.0
        E33          = 15450000000.0
        NU12         = 0.262
        NU13         = 0.264
        NU23         = 0.35
        G12          = 3265000000.0
        G13          = 3495000000.0
        G23          = 3480000000.0

        # Coordinate system
        COORDINATE SYSTEM = sysR
    end parameters for model elastic_orthotropic
end property specification for material glass_uni




  begin solid section hex_section
       formulation = mean_quadrature 
       strain incrementation = midpoint_increment
       #formulation = selective_deviatoric
       #strain incrementation = strongly_objective
  end solid section hex_section

##  Database
  begin finite element model adagio_model
    database name = myBlade_Modified.g
    database type = exodusII
    #decomposition method = rcb

begin parameters for block Adhesive
    material Adhesive
    solid mechanics use model elastic_orthotropic
    section = hex_section
end parameters for block Adhesive


begin parameters for block glass_triax
    material glass_triax
    solid mechanics use model elastic_orthotropic
    section = hex_section
end parameters for block glass_triax


begin parameters for block medium_density_foam
    material medium_density_foam
    solid mechanics use model elastic_orthotropic
    section = hex_section
end parameters for block medium_density_foam


begin parameters for block glass_uni
    material glass_uni
    solid mechanics use model elastic_orthotropic
    section = hex_section
end parameters for block glass_uni



    
  end finite element model adagio_model

  begin adagio procedure Apst_Procedure
    begin time control
        begin time stepping block p1
          start time = 0
          begin parameters for adagio region adagio
            time increment = {tf/150}
          end parameters for adagio region adagio
        end time stepping block p1
        termination time = {tf}
    end time control

###REGION
    begin adagio region adagio
      #begin user variable node_forces_var 
      #  type = node vector
      #end

      use finite element model adagio_model

       ### Results section.


      begin results output results
       database name = myBlade_Modified.e
       database type = exodusii
	   
        start time = {tf}

        # Elastic Orthotropic outputs
         #element variables = unrotated_stress
         #element variables = green_lagrange_strain
         element variables = material_direction_1
         element variables = material_direction_2
         element variables = material_direction_3
         element variables = mat_cauchy_stress_xx
         element variables = mat_cauchy_stress_yy
         element variables = mat_cauchy_stress_zz
         element variables = mat_cauchy_stress_yz
         element variables = mat_cauchy_stress_zx
         element variables = mat_cauchy_stress_xy


         #nodal variables = force_external
         nodal variables = reaction
         nodal variables = displacement
      end results output results
      
# #### initial conditions (material directors)
      #begin initial condition
      #  include all blocks
      #  Initialize Variable Name = node_forces_var
      #  Variable Type = node
      #  Read variable = node_forces
      #end

      Begin Initial Condition
        Include all blocks
        Initialize Variable Name = ANGLE_1
        Variable Type = element
        Read variable = rotation_angle_one
      End

      Begin Initial Condition
        Include all blocks
        Initialize Variable Name = ANGLE_2
        Variable Type = element
        Read variable = rotation_angle_two
      End

      Begin Initial Condition
        Include all blocks
        Initialize Variable Name = ANGLE_3
        Variable Type = element
        Read variable = rotation_angle_three
      End

      Begin Initial Condition
        Include all blocks
        Initialize Variable Name = AXIS_1
        Variable Type = element
        Read variable = rotation_axis_one
      End
	  
      Begin Initial Condition
        Include all blocks
        Initialize Variable Name = AXIS_2
        Variable Type = element
        Read variable = rotation_axis_two
      End
	  
      Begin Initial Condition
        Include all blocks
        Initialize Variable Name = AXIS_3
        Variable Type = element
        Read variable = rotation_axis_three
      End

      
       ################ Kinematic Boundary Conditions ##########
       begin fixed displacement
          node set = station002_ns
          component = x y z
       end
      begin traction
         sideset = oml_ss 
         direction = y
         function = ramp
         scale factor = 1500
      end traction
      begin traction
         sideset = oml_ss 
         direction = x
         function = ramp
         scale factor = -750
      end traction
        ### ------------------###
        ### Solver definition ###
        ### ------------------###
    begin solver
      begin loadstep predictor
         type = secant
      end
      begin cg
         target relative residual = 2.0e-5
         acceptable relative residual = 4.0e-5
         target residual = 1.0e-7
         acceptable residual = 2.0e-7
         maximum iterations       = 3000
         iteration print          = 100
         reference = belytschko 
         preconditioner = probe

        begin full tangent preconditioner
          linear solver = gdsw
          Minimum Smoothing Iterations = 15
          adaptive strategy = update
        end
      end
    end
      begin adaptive time stepping
          cutback factor = 0.5
          growth factor = 1.0
          maximum failure cutbacks = 10
          target iterations = 3000
          iteration window = 200
      end adaptive time stepping
    end adagio region adagio
  end adagio procedure Apst_Procedure

  begin gdsw equation solver gdsw
  end

end Sierra Job