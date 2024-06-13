#include<iostream>
int main(){

	int narray1(3);
	int narray2(3);
	int n_merged = narray1 + narray2;

	int targetArray[6] = {4,5,6,0,0,0};
	int insertArray[3]  = {1,2,3};

	int index_newarray(n_merged-1);
	int index_oriarray(narray1-1);
	int index_insarray(narray2-1);

	while(index_newarray >= 0 || index_oriarray >=0 || index_insarray >= 0){
		if(index_oriarray < 0){
			for(int itmp(index_insarray); itmp >=0; itmp--){
				targetArray[index_newarray] = insertArray[index_insarray];
				index_newarray -= 1;
				index_insarray -= 1;
			}
		}
		else if(index_insarray < 0){
			for(int itmp(index_oriarray); itmp >=0; itmp--){
				targetArray[index_newarray] = targetArray[index_oriarray];
      	index_newarray -= 1;
      	index_oriarray -= 1;
			}
		}
		else{
			if(targetArray[index_oriarray] >= insertArray[index_insarray]){
				targetArray[index_newarray] = targetArray[index_oriarray];
				index_newarray -= 1;
				index_oriarray -= 1;
			}
			else{
				targetArray[index_newarray] = insertArray[index_insarray]; 
				index_newarray -= 1;
				index_insarray -= 1;
			}
		}
	}
	for(unsigned ii(0); ii < n_merged; ii++)std::cout << targetArray[ii] << " ";
	std::cout << std::endl;

	return 1;		
}
