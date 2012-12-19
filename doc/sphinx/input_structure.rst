***************************
Structure of the Input File
*************************** 

Example
=======

The pimc.xml files are xml files. 
The root tag contains six elements. Order doesn't matter to the parser, 
but we always order them in our pimc.xml input files (extended whitespace
indicates omitted details):


.. code-block:: xml

   <?xml version="1.0"?>
   <Simulation>
       <SuperCell     />
       <Species     />
       <Temperature    />
       <Action>
               
       </Action>
       <Estimators>
               
       </Estimators>
       <PIMC    >
               
       </PIMC>
   </Simulation>
    
           
As you can see, the file is grouped by six tags that describe the PIMC 
simulation:

* SuperCell
* Species (this tag is repeated for each species)
* Temperature
* Action (contains all the ActionTags).
* Estimators (contains all the EstimatorTags).
* PIMC (specifies the parallelism and contains the algorithm tree as PIMCTags).
 
All these tags are required.   Notes: The code uses atomic units (Ha, a0), but has some unit conversion capabilities.  

SuperCell
`````````

The SuperCell tag is parsed in SimInfoParser.cc.

Species
```````

The Species tag is parsed in SimInfoParser.cc.

Temperature
```````````

The Temperature tag is parsed in SimInfoParser.cc.

Actions
```````

In path integrals, the action plays the same role as the Hamiltonian plays in Schrödinger's equation. The action section of a pimc.xml file is denoted

              <Actions>
              </Actions>
            
and contains any number of ActionTags. Most simulations contain at least SpringAction for free-particle kinetic action, but in some cases even that may be replaced. Action tags are parsed in ActionParser.cc.

and contains any number of ActionTags. Most simulations contain at least SpringAction for free-particle kinetic action, but in some cases even that may be replaced. Action tags are parsed in ActionParser.cc.

Estimators
``````````

Estimators are the mathematical and algorithmic tools to extract physical information from the path integral. The estimator section is denoted

              <Estimators>
              </Estimators>
          
and may contain any number of EstimatorTags. Estimator tags are parsed in EstimatorParser.cc.

PIMC Commands
`````````````

The PIMC commands are denoted

              <PIMC Commands>
              </PIMC Commands>
          
and may contain any number of sequential or nested PIMCTags. PIMC commands describe the how the paths are sampled, when measurements are performed, and when data is written to disk. PIMC tags are parsed in PIMCParser.cc.
