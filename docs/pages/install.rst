.. _install:

.. contents::
   :depth: 3

=======
Install
=======

To install iGUIDE, simply clone the repository to the desired destination::
  
  git clone https://github.com/cnobles/iGUIDE.git

Then initiate the install using the install script. If you would like the 
installed environment to be named something other than 'iguide', the new conda 
environment name can be provided to the 'install.sh' script as provided below::

  cd path/to/iGUIDE
  bash bin/install.sh

Or::

  cd path/to/iGUIDE
  bash bin/install.sh {env_name}
