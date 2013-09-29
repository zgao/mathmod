#include "GalleryArrangement.h"
#include <stdlib.h>

bool isAllowable(wallList *x){
	
}

void freeArrangement(arrangement *x){
	current = x -> walls;
	while (current != NULL) {
		temp = current;
		freeWall(current -> value);
		current = current -> next;
		free(temp);
	}
}

void freeWall(wall *x){
	current = x;
	while (current != NULL) {
		temp = current;
		current = current -> child;
		free(temp);
	}
}

arrangement* createRandomArrangement(){
	
}
