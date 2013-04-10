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

/**** LOS stands for Linearly Ordered Set  *******/

int SAMEWINDOW=65534;

typedef int word_t;
typedef int index_t;
typedef int keys_t;
typedef int pos_t;

class LOS3 {
    public:
	keys_t id;
	pos_t pos;
	int gap;
        void set(keys_t the_id, pos_t the_pos, int the_gap) {
            id = the_id;
            pos = the_pos;
            gap = the_gap;
        }
        void fprint(FILE *f) {
            fprintf(f," %li:%li:%i",(long)id,(long)pos,gap);
        }
        bool operator<(LOS3 other) {
            return(getid()<other.getid() || (getid()==other.getid() && getpos() < other.getpos()) || (getid()==other.getid() && getpos()==other.getpos() && getgap() > other.getgap()));
        }
        bool operator==(LOS3 other) {
            return(getid()==other.getid() && getpos()==other.getpos() && getgap()==other.getgap());
        }
        void operator=(LOS3 other) {
	    set(other.getid(), other.getpos(), other.getgap());
        }
        bool equiv(LOS3 other) {
            return(getid()==other.getid() && abs(getpos() - other.getpos())<SAMEWINDOW);
        }
	keys_t getid()  { return(id);}
	pos_t  getpos() { return(pos);}
        int    getgap() { return(gap);}
};

#define GAPSHIFT 2
class LOS3A {
    public:
        keys_t id;
        pos_t pos;

        void set(keys_t the_id, pos_t the_pos, int the_gap) {
            id = the_id;
            pos = (the_pos << GAPSHIFT) | the_gap;
        }
        void fprint(FILE *f) {
            fprintf(f,"%li:%li:%i", (long)id, (long)getpos(), getgap());
        }
        bool operator<(LOS3A other) {
            return(getid()<other.getid() || (getid()==other.getid() && getpos() < other.getpos()) || (getid()==other.getid() && getpos()==other.getpos() && getgap() > other.getgap()));
        }
        bool operator==(LOS3A other) {
            return(getid()==other.getid() && getpos()==other.getpos() && getgap()==other.getgap());
        }
        void operator=(LOS3A other) {
            id = other.id;
            pos = other.pos;
        }
	bool equiv(LOS3A other) {
            return(getid()==other.getid() && abs(getpos() - other.getpos())<SAMEWINDOW);
        }
	keys_t getid()  { return(id);}
        pos_t getpos()  { return(pos >> GAPSHIFT);}
        int   getgap()  { return(pos & ((1<<GAPSHIFT)-1));}
};

template <class ordered_set>
class PAIR {
    public:
        ordered_set left;
        ordered_set right;
        void set(ordered_set the_left, ordered_set the_right) {
            left = the_left;
            right = the_right;
        }
        void fprint(FILE *f) {
            fprintf(f,"(");
            left.fprint(f);
            fprintf(f,",");
            right.fprint(f);
            fprintf(f,")");
        }
        void operator=(PAIR other) {
            left = other.left;
            right = other.right;
        }
        bool operator<(PAIR other) {	//LEXICOGRAGHIC ORDER
            return(left<other.left || (left==other.left && right<other.right));
        }
        bool operator==(PAIR other) {	//EQUAL
            return(left==other.left && right==other.right);
        }
        bool equiv(PAIR other) {    	//EQUIVALENCE
            return(left.equiv(other.left) && right.equiv(other.right));
        }
	bool leid(PAIR other) {		//CROSS-ORDER FOR FOLDING
	    return(left.getid()<other.left.getid() || (left.getid()==other.left.getid() && right.getid()<other.right.getid()));
	}

        static void quicksort_less(PAIR arr[], index_t beg, index_t end) {
            index_t i = beg, j = end;
            PAIR tmp;
            PAIR pivot = arr[(beg + end) / 2];
            word_t a;
            while (i <= j) {
                while(arr[i] < pivot) i++;
                while(pivot < arr[j]) j--;
                if (i <= j) {
                    tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
                    i++;
                    j--;
                }
            }
            if (beg < j) quicksort_less(arr, beg, j);
            if (i < end) quicksort_less(arr, i, end);
        }


        static void quicksort_leid(PAIR arr[], word_t word[], index_t beg, index_t end) {
            index_t i = beg, j = end;
            PAIR tmp;
            PAIR pivot = arr[(beg + end) / 2];
            word_t a;
            while (i <= j) {
                while (arr[i].leid(pivot)) i++;
                while (pivot.leid(arr[j])) j--;
                if (i <= j) {
                    tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
                    a = word[i]; word[i] = word[j]; word[j] = a;
                    i++;
                    j--;
                }
            }
            if (beg < j) quicksort_leid(arr, word, beg, j);
            if (i < end) quicksort_leid(arr, word, i, end);
        }

};
