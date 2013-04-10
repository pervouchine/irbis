//	Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
//	This file is a part of the IRBIS package.
//	IRBIS package is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//	
//	IRBIS package is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.


class basepair {
    public:
    	pos_t pos5;
	pos_t pos3;
    	int base5;
	int score;
	pos_t next;

	void set(pos_t the_pos5, pos_t the_pos3, int thebase5){
	    pos5 = the_pos5;
	    pos3 = the_pos3;
	    base5 = thebase5;
	    score = 0;
	    next = -1;
	}

        void fprint(FILE *f) {
            fprintf(f,"(%li:%c-%li:%c %i %i)",(long)pos5, decode1(base5), (long)pos3, decode1(3-base5),score,next);
        }

        bool operator<(basepair other) {
            return(pos5<other.pos5||(pos5==other.pos5 && pos3>other.pos3));
        }

        bool operator==(basepair other) {
            return(pos5==other.pos5 && pos3==other.pos3);
        }

        void operator=(basepair other) {
            pos5 = other.pos5;
            pos3 = other.pos3;
	    base5 = other.base5;
        }

        static void qs(basepair arr[], index_t beg, index_t end) {
            index_t i = beg, j = end;
            basepair tmp;
            basepair pivot = arr[(beg + end) / 2];
            while (i <= j) {
                while (arr[i] < pivot) i++;
                while (pivot < arr[j]) j--;
                if (i <= j) {
                    tmp = arr[i]; 
		    arr[i] = arr[j]; 
		    arr[j] = tmp;
                    i++;
                    j--;
                }
            }
            if (beg < j) qs(arr, beg, j);
            if (i < end) qs(arr, i, end);
        }
	static void compress(basepair arr[], index_t* n) {
	    index_t i,j=0;
	    for(i=1;i<(*n);i++) if(!(arr[i]==arr[j])) arr[++j]=arr[i];
	    (*n) = j+1;
	}
};
