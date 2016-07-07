# Mosquito

This is a re-analysis of some preliminary sequence data from Culex pipiens
form pipiens and form molestus. Individual samples were sequenced at a very
shallow depth. Previous analyses showed a quite low level of mapping success
to the reference genome of Culex quinquefasciatus. Chimeric reads are suspected
to be frequent.

The goal is to compare the two ecotypes.
#2016-06-29
We used SABRE to identificate reads from the same sample(but only for pipiens).
The reads have an adapt which is a barcode, being the one which allows the determination
#2016-06-30
We did 2 tasks: First of all, we used SABRE for form molestus. Secondly, we used a program called PEAR, to merge the reads of C.pipiens
#2016-07-06
We did a unique script with the programs: SABRE and PEAR, to the both forms of the mosquitos 
#2016-07-06b
Now we need to cut off the adapts which are present in our reads. We started working with assembled reads and then continued with unassembled, which probably have more adaptors.
