mathmod
=======

Alright, so I've done some things. 

First, I've added a data structure for an individual in our algorithm
Second, in order to make computation reasonable, we have a fixed arrangement of paintings for a setup of walls. It can be given as follows:

We start by putting one painting on each side in the middle of each portable wall. This part is actually mandated by the problem itself. 

Then, with all of the remaining paintins, we start as close to the camera corners as possible, and add them if they aren't within two meters of a corner of a portable wall.

If we get to 2 meters within of the door on each wall (i.e., we fill up the exterior walls) and we don't have enough paintings, the individual dies. 

This is a pretty big assumption, but it makes a lot of sense. First, there will likely not be very many walls in our setup in the first place, because the cameras should be able to see a lot, or else there will be many dead spots. So cameras will be able to see across the room, or we hope they will be able to. Also, the middle is seen far more often than the sides are, and painints close to the doors can be stolen much more easily anyways, so being directly in the corners is a good thing. 
Plus, since there aren't going to be many walls, the exterior walls will be mostly (if not completely) filled up anyways. So it's not actually that big of an assumption.

Second, this is mostly directed towards MGao, we need to make sure the doors aren't blocked off by external walls. So included in the death check should be that we should have a path from one door to the other. 

Put your code in a .c file with its function prototypes and data structures in a .h file, making sure to use include gaurds. 
(i.e. http://stackoverflow.com/questions/5128664/how-to-split-a-c-program-into-multiple-files) 
