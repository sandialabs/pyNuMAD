// Sierra/SD blank input
   
SOLUTION
   case 'myEigen'
   eigen
   nmodes = 10
END

file
    geometry_file 'myBlade_Modified.g'
end

block medium_density_foam
    material medium_density_foam
end 


block glass_triax
    material glass_triax
end 


block glass_biax
    material glass_biax
end 


block Adhesive
    material Adhesive
end 


block glass_uni
    material glass_uni
end 






outputs
    disp
    stress
    eorient
    material_direction_1
    material_direction_2
    material_direction_3
end
   
material medium_density_foam
orthotropic_prop
    E1          = 142500000.0
    E2          = 142500000.0
    E3          = 142500000.0
    nu12         = 0.3194
    nu23         = 0.3194
    nu13         = 0.3194
    G12          = 54001819.0
    G23          = 54001819.0
    G13          = 54001819.0
    density       = 130.0

end 


material glass_triax
orthotropic_prop
    E1          = 28211400000.0
    E2          = 16238800000.0
    E3          = 15835500000.0
    nu12         = 0.497511
    nu23         = 0.27481
    nu13         = 0.18091
    G12          = 8248220000.0
    G23          = 3491240000.0
    G13          = 3491240000.0
    density       = 1940.0

end 


material glass_biax
orthotropic_prop
    E1          = 11023100000.0
    E2          = 11023100000.0
    E3          = 16047700000.0
    nu12         = 0.688074
    nu23         = 0.117173
    nu13         = 0.117173
    G12          = 13231400000.0
    G23          = 3487480000.0
    G13          = 3487480000.0
    density       = 1940.0

end 


material Adhesive
orthotropic_prop
    E1          = 4560000.0
    E2          = 4560000.0
    E3          = 4560000.0
    nu12         = 0.3
    nu23         = 0.3
    nu13         = 0.3
    G12          = 1450000.0
    G23          = 1450000.0
    G13          = 1450000.0
    density       = 1100.0

end 


material glass_uni
orthotropic_prop
    E1          = 43700000000.0
    E2          = 16500000000.0
    E3          = 15450000000.0
    nu12         = 0.262
    nu23         = 0.35
    nu13         = 0.264
    G12          = 3265000000.0
    G23          = 3480000000.0
    G13          = 3495000000.0
    density       = 1940.0

end 




boundary
   nodeset root
   fixed
   nodeset tip
   z=-10.5
end
