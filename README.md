# CALM
ConservAtion Laws Model

Depedencies
------------------------------

   ROOT
   https://root.cern.ch/drupal/
   cmake

Installation
------------------------------


   mkdir build
   cd build
   cmake ..
   make

   If cmake fails to find ROOT package, you can specify its path via:
   cmake -DROOTSYS='your_path_to_root_installation_e_g_/opt/root'
