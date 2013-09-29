#include "GalleryArrangement.h"
#include <stdlib.h>

bool isAllowable(wallList *x){
	
}

void freeArrangement(arrangement *x){
	wallList *current = x -> walls;
	while (current != NULL) {
		wallList *temp = current;
		freeWall(current -> value);
		current = current -> next;
		free(temp);
	}
}

void freeWall(wall *x){
	wall *current = x;
	while (current != NULL) {
		wall *temp = current;
		current = current -> child;
		free(temp);
	}
}

arrangement* createRandomArrangement(){
	
}
