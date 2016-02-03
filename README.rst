OCOCO - Online Consensus Caller
===============================

Compilation
-----------

.. code-block:: bash

	git clone https://github.com/karel-brinda/ococo
	cd ococo
	cmake .
	make

Compilation for debugging mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Booost detailed logging and asserts.

.. code-block:: bash

	cmake -DDEBUG=ON .
	make

Other options:

* ``VERBOSE_VCF`` - verbose mode for VCF (unchanged bases also reported)
* ``OCOCO32`` - bigger stats (32bits / per position)


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
	Ococo: On-line consensus caller:
	  -i [ --input ] arg      Input SAM/BAM file (- for standard input).
	  -f [ --fasta-ref ] arg  Initial FASTA reference (if not provided, sequence of
	                          N's is considered as the reference).
	  -F [ --fasta-cons ] arg FASTA file with consensus, which is continuously 
	                          updated.
	  -s [ --stats-in ] arg   Input statistics.
	  -S [ --stats-out ] arg  Outputs statistics.
	  -v [ --vcf-cons ] arg   VCF file with updates of consensus.
	  -m [ --mode ] arg       Mode: real-time / batch. [batch]
	  -t [ --strategy ] arg   Strategy for updates: majority / randomized. 
	                          [majority]
	  -a [ --allow-amb ]      Allow updates to ambiguous nucleotides.
	  -q [ --min-MQ ] arg     Minimal mapping quality to increment a counter. [1]
	  -Q [ --min-BQ ] arg     Minimal base quality to increment a counter. [0]
	  -w [ --ref-weight ] arg Initial counter value for nucleotides from the 
	                          reference. [2]
