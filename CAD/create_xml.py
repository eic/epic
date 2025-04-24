#script to create dd4hep xml file
import os
import sys

#initialize some string literals
HEADER = '''
<lccdd>
  <define>
  
    <constant name="SiBarrelMod1_rmin"             value="SiBarrel1_rmin"/> # 269mm
    <constant name="SiBarrelMod2_rmin"             value="SiBarrel2_rmin"/> # 421mm
     <constant name="SiBarrelMod_angle"             value="SiBarrel_angle"/>
    <constant name="SiBarrelMod_dz"                value="SiBarrel_dz"/>

    <constant name="SiBarrelMod1_length"        value="2 * SiBarrelMod1_rmin / tan(SiBarrelMod_angle) - SiBarrel_dz"/>
    <constant name="SiBarrelMod2_length"        value="84*cm"/>

    <constant name="SiBarrelLayer1_length"      value="SiBarrelMod1_length + 1*um"/>
    <constant name="SiBarrelLayer2_length"      value="SiBarrelMod2_length + 1*um"/>
    <constant name="SiBarrelEnvelope_length"    value="SiBarrelLayer2_length + 1*um" />

    <constant name="SiBarrelLayer_thickness"    value="3.0*cm"/>		
    <constant name="SiBarrelLayer1_rmin"        value="SiBarrelMod1_rmin "/>
    <constant name="SiBarrelLayer1_rmax"        value="SiBarrelLayer1_rmin + SiBarrelLayer_thickness"/>
    <constant name="SiBarrelLayer2_rmin"        value="SiBarrelMod2_rmin "/>
    <constant name="SiBarrelLayer2_rmax"        value="SiBarrelLayer2_rmin + SiBarrelLayer_thickness"/>

    <constant name="SiBarrelStaveTilt_angle"     value="0.0*degree"/>
    
  </define>

  <detectors>
    <documentation level="5">
        ### Actual detectors
    </documentation>
    <detector
      id="TrackerBarrel_0_ID"
      name="SagittaSiBarrel"
      type="epic_SiliconBarrelStandardized"
      readout="SiBarrelHits"
      insideTrackingVolume="true">
      <type_flags type="DetType_TRACKER + DetType_BARREL"/>
      <dimensions
        rmin="SiBarrelLayer1_rmin"
        rmax="SiBarrelLayer1_rmax"
        length="SiBarrelLayer1_length" />
      <comment>Silicon Barrel Modules</comment>
      <!-- L3 Stave -->
      <module name="L3Module0" vis="TrackerLayerVis">
        <!--bundle-->
    '''

def string_module_begin(stave_name ,module_number):
    MODULE_BEGIN = f'''<module name="{stave_name}Module{module_number}" vis="TrackerLayerVis">
        <!--bundle-->
    '''
    return MODULE_BEGIN
def string_module_end():
    MODULE_END = '''
        <!--end bundle-->    
      </module> '''
    return MODULE_END


MIDDLE = '''
        <!--end bundle-->    
      </module>
      <comment> Layers composed of many arrayed modules  </comment>
      <layer module="L3Module" id="1" vis="TrackerLayerVis">
        <barrel_envelope
          inner_r="SiBarrelLayer1_rmin-10.0*mm"
          outer_r="SiBarrelLayer1_rmax"
          z_length="SiBarrelLayer1_length + 50*mm" />
        <layer_material surface="inner" binning="binPhi,binZ" bins0="46" bins1="100" />
        <layer_material surface="outer" binning="binPhi,binZ" bins0="46" bins1="100" />
        <comment>
          phi0     : Starting phi of first module.
          phi_tilt : Phi tilt of a module.
          rc       : Radius of the module center.
          nphi     : Number of modules in phi.
          rphi_dr  : The delta radius of every other module.
          z0       : Z position of first module in phi.
          nz       : Number of modules to place in z.
          dr       : Radial displacement parameter, of every other module.
        </comment>
        <rphi_layout phi_tilt="SiBarrelStaveTilt_angle" nphi="46" phi0="0.0" rc="SiBarrelMod1_rmin" dr="7.0 * mm"/>
        <z_layout dr="0.0 * mm" z0="0.0 * mm" nz="1"/>
      </layer>
    </detector>
    <documentation level="5">
        ### Actual detectors
    </documentation>
    <detector
      id="TrackerBarrel_1_ID"
      name="OuterSiBarrel"
      type="epic_SiliconBarrelStandardized"
      readout="SiBarrelHits"
      insideTrackingVolume="true">
      <type_flags type="DetType_TRACKER + DetType_BARREL"/>
      <dimensions
        rmin="SiBarrelLayer2_rmin"
        rmax="SiBarrelLayer2_rmax"
        length="SiBarrelLayer2_length" />
      <comment>Silicon Barrel Modules</comment>
      <comment> L4 Stave </comment>
      <module name="L4Module0" vis="TrackerLayerVis">
        <!--bundle-->
'''
FOOTER = '''
        <!--end bundle-->    
      </module>
      <comment> Layers composed of many arrayed modules  </comment>
      <layer module="L4Module" id="1" vis="TrackerLayerVis">
        <barrel_envelope
          inner_r="SiBarrelLayer2_rmin-10.0*mm"
          outer_r="SiBarrelLayer2_rmax"
          z_length="SiBarrelLayer2_length + 50*mm" />
        <layer_material surface="inner" binning="binPhi,binZ" bins0="128" bins1="100" />
        <layer_material surface="outer" binning="binPhi,binZ" bins0="128" bins1="100" />
        <comment>
          phi0     : Starting phi of first module.
          phi_tilt : Phi tilt of a module.
          rc       : Radius of the module center.
          nphi     : Number of modules in phi.
          rphi_dr  : The delta radius of every other module.
          z0       : Z position of first module in phi.
          nz       : Number of modules to place in z.
          dr       : Radial displacement parameter, of every other module.
        </comment>
        <rphi_layout phi_tilt="SiBarrelStaveTilt_angle" nphi="70" phi0="0.0" rc="SiBarrelMod2_rmin" dr="6.0 * mm"/>
        <z_layout dr="0.0 * mm" z0="0.0 * mm" nz="1"/>
      </layer>
    </detector>
  </detectors>

  <plugins>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="SagittaSiBarrel"/>
      <argument value="layer_pattern: str=SagittaSiBarrel_layer\d"/>
    </plugin>
    <plugin name="DD4hep_ParametersPlugin">
      <argument value="OuterSiBarrel"/>
      <argument value="layer_pattern: str=OuterSiBarrel_layer\d"/>
    </plugin>
  </plugins>

  <readouts>
    <readout name="SiBarrelHits">
      <segmentation type="CartesianGridXY" grid_size_x="0.020*mm" grid_size_y="0.020*mm" />
      <comment>bitfieldsizes for indexing sensors and modules, max values are 2 exp bitfieldsize</comment>
      <id>system:8,layer:1,module:11,sensor:9,x:32:-16,y:-16</id>
    </readout>
  </readouts>

</lccdd>

'''

MAX_COMPONENT_NUMBER = 80 #less then 2 to the 9 which is what we allocated for sensors

#note that units are in milimeters
default_thickness_sensitive = 0.09996 #in milimeters
def module_component(component_path, component_name, my_dict, file_count, stave_name):
    COMPONENT_HEADER = f'''
        <module_component name="{component_name}" '''

    sens = False
    #print the attributes of the dictionary
    COMPONENT_BODY = ""      
    for key, value in my_dict.items():
        COMPONENT_BODY += f"\n                          {key}=\"{value}\""
                          #material="{component_material}"
                          #sensitive="{component_offset}"
                          #offset="{component_is_sensitive} * mm"
                          #thickness="{component_thickness} * mm"
        if key == "sensitive" and value == "true":
          sens = True
    COMPONENT_FOOTER = f'''
                          vis="TrackerLayerVis"
                          file="{component_path}" />
    '''
    COMPONENT = COMPONENT_HEADER + COMPONENT_BODY + COMPONENT_FOOTER
    
    #note: do module nesting only for sensitive components (note that this approach might introduce some bugs)
    remainder = file_count % MAX_COMPONENT_NUMBER
    if remainder + 1 == MAX_COMPONENT_NUMBER and sens:
      module_number = int(file_count/MAX_COMPONENT_NUMBER)
      COMPONENT = string_module_end() + "\n    " + string_module_begin(stave_name ,module_number + 1) + COMPONENT
      return COMPONENT
    else:
      return COMPONENT
    
file_count_global = 0
def scan_and_place(folder_path, search_string, xml_file, my_dict, stave_name):
    # scan folder for gdml name with string and place certain attributes accordingly
    # i.e the material etc. of a component is determined by the name of the .gdml file (very lazy way to do it I know)
    #add the counter HERE
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".gdml") and search_string in file:
                # Get the relative path and add to the list
                full_path = os.path.join(root, file)
                global file_count_global
                xml_file.write(module_component(full_path, file, my_dict, file_count_global, stave_name))     
                file_count_global += 1



#dictionaries to contain the filename patterns
dict_list = [
    {"matching_name": "Active", "material": "Silicon", "sensitive": "true", "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "ALICE", "material": "Silicon", "sensitive": "true", "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "Biasing", "material": "Silicon"      , "sensitive": "false" , "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "DataBackbone", "material": "Silicon" , "sensitive": "false" , "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "Pads", "material": "Silicon"         , "sensitive": "false" , "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "PowerSwitches", "material": "Silicon", "sensitive": "false" , "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "Readout", "material": "Silicon"      , "sensitive": "false" , "thickness": f"{default_thickness_sensitive} * mm"},
    {"matching_name": "FPC", "material": "Kapton"           , "sensitive": "false"},
    {"matching_name": "kapton", "material": "Kapton"        , "sensitive": "false"},
    {"matching_name": "Kapton", "material": "Kapton"        , "sensitive": "false"},
    {"matching_name": "K9", "material": "K9"                , "sensitive": "false"},
    {"matching_name": "Carbon", "material": "CarbonFiber"   , "sensitive": "false"},
    {"matching_name": "Ultem", "material": "Ultem"          , "sensitive": "false"},
] #extend as required
dict_attr_list = ["material", "sensitive", "thickness", "offset"]
  

# Check if the folder path is provided as an argument
if len(sys.argv) < 3:
    print("Please provide the L3 and L4 gdml folder paths as the first and second command line arguments respectively.")
    sys.exit(1)

L3_folder_path = sys.argv[1]
L4_folder_path = sys.argv[2]
file_path = "silicon_barrel.xml"
# Check if the file exists
if os.path.exists(file_path):
    raise FileExistsError(f"Error: The file '{file_path}' already exists.")
    sys.exit(1)

with open(file_path, "w") as xml_file:
    xml_file.write(HEADER)

    file_count_global = 0 #reset the global file count in each new stave
    for key_dict in dict_list:   
        scan_and_place(L3_folder_path, key_dict["matching_name"], xml_file, key_dict, "L3")
    xml_file.write(MIDDLE)
    file_count_global = 0
    for key_dict in dict_list:
        scan_and_place(L4_folder_path, key_dict["matching_name"], xml_file, key_dict, "L4")
    xml_file.write(FOOTER)

# Print the resulting list (optional)
#print("Found .gdml files:")
#for file in gdml_files:
#    print(file)
