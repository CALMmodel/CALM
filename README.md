# CALM
ConservAtion Laws Model

Depedencies
------------------------------

   cmake
   
   ROOT (http://root.cern.ch/)

Installation
------------------------------

      mkdir build
      cd build
      cmake ..
      make
   
   
   If cmake fails to find ROOT package, you can specify its path via:
   
      cmake -DROOTSYS='your_path_to_root_installation_e_g_/opt/root'
