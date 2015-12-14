OCOCO - Online Consensus Caller
===============================

Compilation
-----------

.. code-block:: bash

	git clone https://github.com/karel-brinda/ococo
	cd ococo
	cmake .
	make


Installation
------------

.. code-block:: bash
	
	make install


Tests
-----

.. code-block:: bash

	ctest --output-on-failure


Parameters
----------

.. code-block::

	./ococo 
	Ococo: On-line consensus calling.:
	  -i [ --input ] arg           Input SAM/BAM file (- for standard input).
	  -f [ --init-fasta ] arg      Initial FASTA reference.
	  -c [ --consensus-fasta ] arg Consensus (up-to-date FASTA reference).
	  -s [ --stats ] arg           File with up-to-date statistics.
	  -m [ --min-map-qual ] arg    Minimal mapping quality. [1]
	  -b [ --min-base-qual ] arg   Minimal base quality. [0]

