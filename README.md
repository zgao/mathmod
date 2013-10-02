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

Also, for the death check, I reread the problem spec and Parallel walls are not allowed to be within 5 (or 4.99999 because floats) meters of eachother.
So we simply check to see for every corner if there is a corner less than 5 meters from it. That makes it easy. So much eaier than intersection checking. It also accounts for if a wall is close but not remotely "parallel" (imagine two circles of radius 5m around each endpoint of a wall. 

I added some of the stuff for deallocating individuals, someone needs to write the death check and the generation of the random individual. Maybe me. 

Put your code in a .c file with its function prototypes and data structures in a .h file, making sure to use include gaurds. 
(i.e. http://stackoverflow.com/questions/5128664/how-to-split-a-c-program-into-multiple-files) 


Also, we don't consider tree-like walls. So each wall is only connected to at most one other on each end. This makes computation a LOT simpler, and also, having a 3 way wall necessarily creates deadspots (unless the wall is pointing exactly at the camera, but in reality that is actually just 2 dead spots, because you can't see the front of either painting). 

On another note, the data structure is really hackey. It allows for one point walls, and the last element has an angle but no wall attached to it. Think of the wall structus as posts that reside on the corners, and have one direction for walls to attach to it. If there is no child then it's a terminal post. There is a method that returns the points of the wall, it's called wallPoints in GalleryArrangement.c
-Sam 
