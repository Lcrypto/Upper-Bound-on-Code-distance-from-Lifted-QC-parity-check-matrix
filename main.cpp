#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <NTL/GF2X.h>
#include <stdio.h>
#include <iomanip>
#include <cassert>
#include "..\src\include\mpi.h"
#include <algorithm>

NTL_CLIENT
using namespace std;

long nchoosek(int n, int k) {
	double res = 1;
	double d_n = n;
	double d_k = k;
	for (int i = 1; i <= k; ++i) {
		res *= (d_n - (d_k - i)) / i;
	}
	return res;
}


GF2X divider = GF2X();
void print(const GF2X& poly) {
	int length = deg(poly);
	if (length == -1) {
		cout << "Empty poly" << endl;
		return;
	}
	for (int i = 0; i < length + 1; ++i) {
		cout << coeff(poly, i);
	}
	cout << endl;
}

GF2X poly_mult(const GF2X& poly1, const GF2X& poly2) {
	GF2X tmp = poly1 * poly2;
	GF2X r;
	rem(r, tmp, divider);
	r.normalize();
	return r;
}

long poly_weight(const GF2X& poly) {
	return weight(poly);
}

class Node {
	int start;

	int children_count;
	int value;
	Node** children;
public:
	Node() {
		start = -1;
		children_count = -1;
		value = -1;
		children = 0;
	}

	void init(int start, int step, int rows, int columns) {

		int max_count = columns - rows + step + 1;
		if (max_count <= start) {
			return;
		}
		if (step == rows) {
			return;
		}
		children_count = max_count - start;
		children = new Node*[children_count];
		this->start = start;
		for (int i = 0; i < children_count; ++i) {
			children[i] = new Node();
		}
	}
	int get_count() {
		return children_count;
	}


	Node* get_item(int ind) {
		return children[ind - start];
	}

	void set_value(int value) {
		this->value = value;
	}

	int get_value() {
		return value;
	}

	~Node() {
		//cout << children << endl;
		if (children == 0)
			return;
		//cout << children_count << endl;
		for (int i = 0; i < children_count; ++i) {
			delete children[i];
		}
		delete[] children;
	}
};

class PermanentStorage {
	Node* root;
	int rows, columns;

	void init(Node* child, int start, int step, int rows, int columns) {
		//child = new Node();
		child->init(start, step, rows, columns);
		int max_count = child->get_count() + start;
		for (int i = start; i < max_count; ++i) {
			init(child->get_item(i), i + 1, step + 1, rows, columns);
		}
	}


public:
	PermanentStorage() {
		root = 0;
	}
	PermanentStorage(int rows, int columns) {
		root = new Node();
		init(root, 0, 0, rows, columns);
	}

	void set_value(const vector<int>& combination, int value) {
		Node* cur = root;
		for (int i = 0; i < combination.size(); ++i) {
			cur = cur->get_item(combination[i]);
		}
		cur->set_value(value);
	}

	void set_value(const vector<int>& combination, int erased_position, int value) {
		Node* cur = root;
		for (int i = 0; i < combination.size(); ++i) {
			if (i == erased_position) {
				continue;
			}
			cur = cur->get_item(combination[i]);
		}
		cur->set_value(value);
	}

	int get_value(const vector<int>& combination) {
		Node* cur = root;
		for (int i = 0; i < combination.size(); ++i) {
			cur = cur->get_item(combination[i]);
		}
		return cur->get_value();
	}


	int get_value(const vector<int>& combination, int erased_position) {
		Node* cur = root;
		for (int i = 0; i < combination.size(); ++i) {
			if (i == erased_position) {
				continue;
			}
			cur = cur->get_item(combination[i]);
		}
		return cur->get_value();
	}

	~PermanentStorage() {
		if (root == 0)
			return;
		delete root;
	}

};

int lim_value = 1000050000;
int global_min = lim_value;
int combination[256];
int rows, columns, circ_size;
int non_zero_elements[128][128];
int non_zero_elements_in_row[128];
float sum_time;
int permutation_flags[128];

GF2X get_permanent_recursive(GF2X& current_value, int step, const vector<vector<GF2X> >& mtr) {
	GF2X res;
	res.SetLength(circ_size);
	if (step >= rows) {
		return current_value;
	}
	for (int i = 0; i < non_zero_elements_in_row[step]; ++i) {
		if (permutation_flags[non_zero_elements[step][i]]) {
			continue;
		}
		permutation_flags[non_zero_elements[step][i]] = 1;

		GF2X t = poly_mult(current_value, mtr[non_zero_elements[step][i]][step]);
		res += get_permanent_recursive(t, step + 1, mtr);
		permutation_flags[non_zero_elements[step][i]] = 0;
	}
	return res;
}

int get_permanent(const vector<int>& mask, const vector<vector<GF2X> >& mtr, const vector<vector<int> >& proto_mtr) {
	int ind[128];
	for (int i = 0; i < rows; ++i) {
		ind[i] = 0;
	}
	for (int c = 0; c < rows; ++c) {
		for (int r = 0; r < rows; ++r) {

			if (proto_mtr[mask[c]][r]) {
				non_zero_elements[r][ind[r]++] = mask[c];
			}
		}

	}

	for (int i = 0; i < rows; ++i) {
		if (ind[i] == 0)
			return 0;
		non_zero_elements_in_row[i] = ind[i];
	}

	clock_t t1 = clock();
	GF2X initial_poly = GF2X();
	initial_poly.SetLength(1);
	initial_poly[0] = 1;
	GF2X res = get_permanent_recursive(initial_poly, 0, mtr);
	clock_t t2 = clock();
	float diff = ((float)t2 - (float)t1);
	sum_time += diff;
	return weight(res);
}

int solve_my(const vector<int>& mask, const vector<vector<GF2X> >& mtr, vector<vector<int> >& proto_mtr, PermanentStorage& storage) {
	long res = 0;
	for (int i = 0; i < mask.size(); ++i) {
		int stored_value = storage.get_value(mask, i);
		if (stored_value != -1) {
			res += stored_value;
		}
		else {
			stored_value = get_permanent(mask, mtr, proto_mtr);
			storage.set_value(mask, i, stored_value);
			res += stored_value;
		}
		if (res >= global_min) {
			return 0;
		}
	}
	return res;
}


void check_all_combinations(int step, int beg, const vector<vector<GF2X> >& mtr, vector<vector<int> >& proto_mtr, PermanentStorage& storage) {
	if (step == rows + 1) {
		vector<int> mask = vector<int>(combination, combination + rows + 1);
		long res = solve_my(mask, mtr, proto_mtr, storage);
		if (res >= global_min)
			return;

		if (res != 0) {
			//cout<<"new min: "<<res<<endl;
			global_min = res;
		}
		return;
	}
	for (int i = beg; i < columns; ++i) {
		combination[step] = i;
		check_all_combinations(step + 1, i + 1, mtr, proto_mtr, storage);
	}
}




void run(const vector<vector<GF2X> >& mtr, vector<vector<int> >& proto_mtr, PermanentStorage& storage) {
	check_all_combinations(0, 0, mtr, proto_mtr, storage);
}

bool nextCombination(vector<int>& a, int n) {
	int k = (int)a.size();
	for (int i = k - 1; i >= 0; --i)
		if (a[i] < n - k + i + 1) {
			++a[i];
			for (int j = i + 1; j < k; ++j)
				a[j] = a[j - 1] + 1;
			return true;
		}
	return false;
}

int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	// read matrix, initialize global variables
	vector<pair<int, string> > results;
	string filelist_filename = "data_for_check1.txt";//"data_for_check1.txt";
	string output_filename = "data_for_check1_out.txt";//"data_for_check1_out.txt";

	for (int i = 1; i < argc; ++i) {
		if (string(argv[i]) == "-inputFile") {
			filelist_filename = argv[i + 1];
			++i;
			continue;
		}
		if (string(argv[i]) == "-outputFile") {
			output_filename = argv[i + 1];
			++i;
			continue;
		}

	}
	global_min = lim_value;
	if (filelist_filename == "") {
		cerr << "wrong input\n";
		//cerr << "Usage: DFreeUpperBoundOnProtograph.exe -inputFile INPUT.TXT -outputFile OUTPUT.TXT\n";
		return 1;
	}
	clock_t t1, t2;
	ofstream output_file;
	if (ProcRank == 0) {
		output_file.open(output_filename, ios::app);
		t1 = clock();
	}
	std::ifstream filelist_input(filelist_filename);
	bool is_odd = false;
	for (string filename; getline(filelist_input, filename);) {
		is_odd ^= 1;
		if (ProcRank == 0) {
			cout << filename << endl;
		}
		std::ifstream input(filename);

		string fstline;
		getline(input, fstline);
		istringstream iss(fstline);
		if (!(iss >> columns >> rows >> circ_size)) {
			if (ProcRank == 0) {
				cout << "Error occured while reading input file" << endl;
			}
		}

		vector<vector<GF2X> > checkmatrix_poly(columns, vector<GF2X>(rows));
		for (int c = 0; c < columns; ++c) {
			for (int r = 0; r < rows; ++r) {
				checkmatrix_poly[c][r].SetLength(circ_size);
			}
		}

		int r = 0;
		for (int r = 0; r < rows; ++r)
		{
			string line;
			getline(input, line);
			istringstream iss(line);
			string matrix_entry;
			for (int c = 0; c < columns; ++c) {
				iss >> matrix_entry;
				if (matrix_entry.find('&') != -1) {
					istringstream iss_entry(matrix_entry);
					string token;

					while (getline(iss_entry, token, '&')) {
						int x = stoi(token);
						checkmatrix_poly[c][r][x] = 1;
					}
				}
				else {
					int x = stoi(matrix_entry);
					if (x == -1) {
						continue;
					}
					checkmatrix_poly[c][r][x] = 1;

				}
				checkmatrix_poly[c][r].normalize();
			}
		}
		vector<vector<int> > protomatrix(columns, vector<int>(rows));
		for (int c = 0; c < columns; ++c) {
			for (int r = 0; r < rows; ++r) {
				protomatrix[c][r] = weight(checkmatrix_poly[c][r]);
			}
		}

		/*for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < columns; ++c) {
		cout << protomatrix[c][r];
		}
		cout << endl;
		}*/
		divider.SetLength(circ_size + 1);
		divider[0] = 1;
		divider[circ_size] = 1;

		//permanents calculation
		int combinations_count = nchoosek(columns, rows);
		int** all_combinations = new int*[combinations_count];
		vector<int> cur_combination(rows);
		for (int i = 0; i < rows; ++i) {
			cur_combination[i] = i;
		}
		int ind = 0;
		do {
			all_combinations[ind] = new int[rows];
			for (int i = 0; i < rows; ++i) {
				all_combinations[ind][i] = cur_combination[i];
			}
			++ind;
		} while (nextCombination(cur_combination, columns - 1));
		if (ind != combinations_count) {
			cerr << "N choose k return wrong value (in fact, it was predictable)" << endl;
			return -1;
		}

		int fake_proc_rank = ProcRank;
		if (!is_odd && ProcRank) {
			fake_proc_rank = ProcNum - ProcRank;
		}
		int perm_count_per_proc = combinations_count / ProcNum;
		int start_perm_ind = perm_count_per_proc*fake_proc_rank;
		int end_perm_ind = perm_count_per_proc*(fake_proc_rank + 1);
		if (fake_proc_rank == ProcNum - 1) {
			end_perm_ind = combinations_count;
		}

		int* permanents = new int[end_perm_ind - start_perm_ind];
		if (ProcRank == 0) {
			cout << "Each proc computes " << end_perm_ind << " permanents" << endl;
		}
		for (int i = start_perm_ind; i < end_perm_ind; ++i) {
			
			vector<int> mask(all_combinations[i], all_combinations[i] + rows);
			permanents[i - start_perm_ind] = get_permanent(mask, checkmatrix_poly, protomatrix);
		}

		
		if (ProcRank == 0) {
			int* permanents_from_all = new int[combinations_count];
			for (int i = 0; i < end_perm_ind; ++i) {
				permanents_from_all[i] = permanents[i];
			}
			MPI_Status status;
			for (int p = 1; p < ProcNum; ++p) {
				fake_proc_rank = p;
				if (!is_odd) {
					fake_proc_rank = ProcNum - p;
				}
				int start_perm_ind = perm_count_per_proc*fake_proc_rank;
				int end_perm_ind = perm_count_per_proc*(fake_proc_rank + 1);
				if (fake_proc_rank == ProcNum - 1) {
					end_perm_ind = combinations_count;
				}
				MPI_Recv(permanents_from_all + start_perm_ind, end_perm_ind - start_perm_ind, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
				//cout << "Permanents values received form proc " << p << endl;
			}

			PermanentStorage storage = PermanentStorage(rows, columns);
			for (int i = 0; i < combinations_count; ++i) {
				vector<int> mask(all_combinations[i], all_combinations[i] + rows);
				storage.set_value(mask, permanents_from_all[i]);
			}
			run(checkmatrix_poly, protomatrix, storage);
			cout << "Upper bound: " << global_min << endl;
			output_file << global_min << " " << filename << endl;
			global_min = lim_value;
			delete[] permanents_from_all;
		}
		else {
			MPI_Send(permanents, end_perm_ind - start_perm_ind, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}

		delete[] permanents;
		for (int i = 0; i < combinations_count; ++i) {
			delete[] all_combinations[i];
		}
		delete[] all_combinations;
	}



	if (ProcRank == 0) {
		output_file.close();
		t2 = clock();
		float diff = ((float)t2 - (float)t1);
		cout << diff / CLOCKS_PER_SEC << endl;
	}


	MPI_Finalize();
	return 0;
}


