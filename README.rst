Kidscpp
=======

Kidscpp processes the kinetic inductance detector (KID) signals for
`TolTEC <http://toltec.astro.umass.edu>`_.

It is developed as part of TolTECA, the TolTEC data analysis software suite.

While kidscpp is developed targeting TolTEC KIDs and the ROACH-2 readout
system, it can also be adapted to other instruments that shares similar
architectural properties.


Features
--------

Kidscpp can:

* read raw KIDs data files generated from the data acquisition system.

* find KID resonance peaks in multiplexed frequency sweep data.

* fit KID resonance model in multiplexed frequency sweep data.

* solve raw KID signal stream (I, Q) using the best fit KIDs model to get
  stream (r, x), where x is the optical detuning proportional to the
  change of detector power loading, and r is the noise term.

* write process result to files.


Kidscpp can be used either as a standalone program (the :code:`kids`
executable) or as a library to be integrated in a larger data analysis pipeline
(the :code:`libkids` lib).

Build
-----

Kidscpp shares the same build requirements and procedures as
`citlali <https://github.com/toltec-astro/citlali/tree/v0.1.x>`_.

Please refer to the build instructions therein.


Usage
-----

Once successfully built, the created executables will be available in
:code:`build/bin`.

To check the version of the program:

.. code-block::

    // In the build directory:
    $ ./bin/kids --version

To show the help screen of the commandline interface:

.. code-block::

    // In the build directory:
    $ ./bin/kids --help

Please see the `API documentation
<https://toltec-astro.github.io/kidscpp>`_ for details.


License
-------

3-Clause BSD.
