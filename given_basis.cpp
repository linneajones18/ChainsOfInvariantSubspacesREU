#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <cmath>

using namespace std;
using namespace Eigen;

/*
    Instructions for running program:
    1. open a terminal and navigate to the directory this file is in
    2. ensure you've downloaded eigen-3.4.0 from the internet to this directory
    3. enter the following command into the terminal to compile:
        g++ -std=c++17 given_basis.cpp -I eigen-3.4.0/
    4. enter the following command into the terminal to run:
        ./a.out
*/



// -------------------------------------  HELPER FUNCTIONS   ---------------------------------------

// loading bar
// taken from chatgpt
void printLoadingBar(int progress, int total, int barWidth = 50) {
    float fraction = static_cast<float>(progress) / total;
    int filledWidth = static_cast<int>(barWidth * fraction);

    // Build the loading bar string
    string bar = "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < filledWidth) {
            bar += "=";
        } else {
            bar += " ";
        }
    }
    bar += "] " + to_string(int(fraction * 100)) + "%";

    // Print the loading bar
    cout << "\r" << bar;  // Print the loading bar with a carriage return
    cout.flush();         // Ensure the output is displayed immediately
}

// converts an int to a string of its binary form
// taken from internet
string decimalToBinary(int num) {
    if (num == 0) {
        return "0";
    }

    string binaryStr;
    while (num > 0) {
        int remainder = num % 2;
        binaryStr += to_string(remainder);
        num /= 2;
    }

    // The binary string is generated in reverse order
    reverse(binaryStr.begin(), binaryStr.end());

    return binaryStr;
}

// gets a set of vectors from user
vector<VectorXd> getBasis()
{
    vector<VectorXd> basis; //stores our basis
    int length; //length of all the vectors in the basis

    cout << "Would you like to (1) manually enter your vectors or (2) import them from a file?" << endl;
    string choice;
    cin >> choice;
    while((choice != "1") && (choice != "2"))
    {
        cout << "Please enter 1 or 2" << endl;
        cin >> choice;
    }


    if(choice == "1")
    {
        string input = ""; //stores our user input
        bool invalid; //stores the validity status of the vector the user entered
        vector<bool> match;

        // take user input as a string and put each element in the string in a vector
        cout << "Please enter the vectors in the basis, enter 'DONE' when finished:" << endl;

        // takes the first vector input over and over until the user enters a sequence that only consists of 0s and 1s
        do {
            invalid = false;
            input = "";
            cin >> input;
            for(int i = 0; i < input.length(); i++)
            {
                if((input[i] != '0') && (input[i] != '1'))
                {
                    invalid = true;
                    break;
                }
            }

            if(invalid)
            {
                cout << "The input you entered is not exclusively composed of 0s and 1s, please try again:" << endl;
            }

        } while(invalid);

        // add first vector to basis (outside of the next while loop because this decides the lenth of all the other vectors)
        length = input.length();
        VectorXd v(length);
        for(int i = 0; i < length; i++)
        {
            v[i] = input [i];
            match.push_back(false);
        }
        basis.push_back(v);

        // now take the user input for the next vector
        cout << "Vector added successfully, please enter the next vector or enter 'DONE' if finished:" << endl;
        input = "";
        cin >> input;

        // repeat this process of adding vectors until the user decides they're finished
        while((input != "DONE") && (input != "done"))
        {
            invalid = false;
            if(input.length() != length)
            {
                cout << "The vector you added does not match the length of the previously added vector(s). Please try again or enter 'DONE' if you're finished." << endl;
                input = "";
                cin >> input;
                continue;
            }

            // check that the new vector matches the size of the previous ones
            for(int i = 0; i < length; i++)
            {
                if((input[i] != '0') && (input[i] != '1'))
                {
                    invalid = true;
                    break;
                }

                v[i] = input[i];
            }
            if(invalid)
            {
                invalid = false;
                cout << "The input you entered is not exclusively composed of 0s and 1s, please try again or enter 'DONE' if you're finished:" << endl;
                input = "";
                cin >> input;
                continue;
            }

            // make sure vector is not already in the database:
            for(int i = 0; i < basis.size(); i++)
            {
                // reset our match var to assume true until found to be false
                invalid = true;

                // reset the match vector
                for(int j = 0; j < length; j++)
                {
                    match.at(j) = false;
                }

                // this loop goes through the ith sequence in the basis and creates a bool array that tells us if jth character in the two vectors match
                for(int j = 0; j < length; j++)
                {
                    if(basis.at(i)[j] == input[j])
                    {
                        match.at(j) = true;
                    }
                }

                // this loop goes through the array we just made to see if all values are true - if they are - this vector already exists in the basis and we exit the while loop
                for(int j = 0; j < length; j++)
                {
                    if(match.at(j) == false)
                    {
                        invalid = false;
                        break;
                    }
                }
                if(invalid)
                {
                    break;
                }
            }
            if(invalid)
            {
                cout << "You have already entered this vector, please enter another or enter 'DONE' if you're finished:" << endl;
                input = "";
                cin >> input;
                continue;
            }

            // if everything about the new vector checks out - add it to the basis
            basis.push_back(v);
            cout << "Vector added successfully, please enter the next vector or enter 'DONE' if finished:" << endl;
            input = "";
            cin >> input;

        }

        // erase match vector for memory purposes
        match.clear();
    }
    
    else
    {
        ifstream file;
        cout << "What is the name of the file you'd like to use?" << endl;
        string file_name;
        cin >> file_name;
        file.open(file_name);
        if(file.fail())
        {
            cout << "Failed to find a file named: " << file_name << endl << "Returning an empty basis." << endl;
            return basis;
        }

        string cur_line;
        while(!file.eof())
        {
            getline(file, cur_line);
            if(!cur_line.empty())
            {
                length = cur_line.length();
                VectorXd v(length);
                for(int i = 0; i < length; i++)
                {
                    v[i] = cur_line[i];
                }
                basis.push_back(v);
            }
        }

        file.close();
    }

    //fix bug that adds 48 to each entry:
    for(int i = 0; i < basis.size(); i++)
    {
        for(int j = 0; j < length; j++)
        {
            if((basis.at(i)[j] != 1) && (basis.at(i)[j] != 0))
            {
                basis.at(i)[j] -= 48;
            }
        }
    }

    return basis;
}

// reduces any sized matrix to echelon form
MatrixXd rowReduce(MatrixXd mat) 
{
    const int rows = mat.rows();
    const int cols = mat.cols();
    bool pivot_exists = false;
    int starting_row = 0;

    /** Algorithm:
     * find a row that contains a one in the current column
     * swap ith row, where i is the current column
     * add this row to every row below it that contains a one in the ith column
    */

   for(int i = 0; i < cols; i++)
   {
        //if there was no pivot found in the last column, the starting row needs to stay the same, otherwise we move on to the next row
        if(pivot_exists)
        {
            starting_row++;
        }
        pivot_exists = false;

        //for matrices that have more columns than rows, this will prevent it from trying to find more pivot columns when it's already gone through all the rows
        if(starting_row == rows - 1)
        {
            break;
        }

        //this loop searches for a pivot row and swaps it to bring it to the top
        for(int j = starting_row; j < rows; j++)
        {
            if(mat(j, i) == 1)
            {
                mat.row(starting_row).swap(mat.row(j));
                pivot_exists = true;
                break;
            }
        }

        //if there was not a pivot row found for this column, move on to the next column (next iteration of the loop)
        if(!pivot_exists)
        {
            starting_row--;
            continue;
        }

        // goes through all the rows after the new pivot row we just swapped
        for(int j = rows - 1; j > starting_row; j--)
        {
            // if current row has a pivot in col i, add the pivot row to this current row
            if(mat(j, i) == 1)
            {
                // goes through each col of the row
                for(int k = 0; k < cols; k++)
                {
                    if(((mat(starting_row, k) == 0) && (mat(j, k) == 0)) || (mat(starting_row, k) == 1) && (mat(j, k) == 1))
                    {
                        mat(j, k) = 0;
                    }
                    else
                    {
                        mat(j, k) = 1;
                    }
                }
            }
        }
   }

   return mat;
}

// returns the rank of the passed matrix
int getRank(MatrixXd mat)
{
    int rank = 0;
    int cols = mat.cols();
    int rows = mat.rows();
    mat = rowReduce(mat);

    // counts the # of pivot columns (# of nonzero rows)
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            if(mat(i, j) != 0)
            {
                rank++;
                break;
            }
        }
    }

    return rank;
}

// returns true if the given vector of vectors form a basis - false otherwise
bool isValidBasis(vector<VectorXd> basis)
{
    //an easy check that there aren't too many vectors to be a basis, to shorten runtime possibly:
    int length = basis.at(0).rows();
    int basis_size = basis.size();

    if(basis_size > length)
    {
        return false;
    }

    // turn the basis into a matrix
    MatrixXd mat(length, basis_size);

    for(int i = 0; i < basis_size; i++)
    {
        mat.col(i) = basis.at(i);
    }

    int rank = getRank(mat);
    if(rank == mat.cols())
    {
        return true;
    }

    return false;
}

// returns true if the given vector of vectors form an invariant basis - false otherwise
bool isInvariantBasis(vector<VectorXd> basis)
{
    if(!isValidBasis(basis))
    {
        return false;
    }

    int basis_size = basis.size();
    int length = basis.at(0).rows();
    MatrixXd mat(length, 2 * basis_size);
    vector<VectorXd> cycled_vecs;
    VectorXd v(length);
    int last;

    //cycle the vectors
    for(int i = 0; i < basis_size; i++)
    {
        v = basis.at(i);
        last = v(length - 1);
        for(int j = length - 1; j > 0; j--)
        {
            v(j) = v(j - 1);
        }
        v(0) = last;
        cycled_vecs.push_back(v);
    }

    //add the original basis to the matrix
    for(int i = 0; i < basis_size; i++)
    {
        mat.col(i) = basis.at(i);
    }

    //add the cycled basis to the matrix
    for(int i = 0; i < basis_size; i++)
    {
        mat.col(i + basis_size) = cycled_vecs.at(i);
    }

    int rank = getRank(mat);
    cycled_vecs.clear();

    //we already know we have a basis, for it to be invariant, the rank of the matrix must be the same with the cycled vectors as without
    if(rank == basis_size)
    {
        return true;
    }


    return false;
}

// checks if two bases represent different chains of subspaces
bool areDiffSpaces(vector<VectorXd> v1, vector<VectorXd> v2)
{
    if(v1 == v2)
    {
        return false;
    }

    MatrixXd mat1(4, 2);
    mat1 << v1.at(0), v2.at(0);
    if(getRank(mat1) != 1)
    {
        return true;
    }

    if(v1.size() == 1)
    {
        return false;
    }

    MatrixXd mat2(4, 4);
    mat2 << mat1, v1.at(1), v2.at(1);
    if(getRank(mat2) != 2)
    {
        return true;
    }

    if(v1.size() == 2)
    {
        return false;
    }

    MatrixXd mat3(4, 6);
    mat3 << mat2, v1.at(2), v2.at(2);
    if(getRank(mat3) != 3)
    {
        return true;
    }

    return false;
}




// -----------------------  FUNCTIONS FOR SOLVING ORIGINAL CODE JAM PROBLEM   -------------------------

// creates a basis for an inputted dimension and puts it in a file
void createBasis(int dim)
{
    //first check if its a power of 2:
    if(!((dim != 0) && ((dim & (dim - 1)) == 0)))
    {
        cout << "The dimension entered is not a power of 2" << endl;
        return;
    }

    ifstream file;
    string file_name = to_string(dim) + "bit.txt";
    file.open(file_name);

    if(!file.fail())
    {
        cout << "File already exists" << endl;
        return;
    }

    /** algo:
     * see if we have a file for the dimension down - if not, we're going to make one
     * use the previous dimension file to make this new one
    */

   //find the next availale file below
   int dim_below = dim;
   do
   {
        dim_below = dim_below / 2;
        file_name = to_string(dim_below) + "bit.txt";
        file.open(file_name);
   } while ((file.fail()) && (dim_below > 2));

        if(dim_below == 2)
        {
            ofstream first_file;
            first_file.open("4bit.txt");
            first_file << "1111" << "\n";
            first_file << "0101" << "\n";
            first_file << "0011" << "\n";
            first_file << "0001" << "\n";
            first_file.close();
        }

        if(dim == 4)
        {
            return;
        }
   
    //open the lowest file and create as many files as necessary to get to the input one
    while(dim_below < dim)
    {
        //load the contents of the file below into a basis vec
        ifstream file_below;
        file_below.open(to_string(dim_below) + "bit.txt");
        int length;
        vector<VectorXd> basis;
        string cur_line;

        while(!file_below.eof())
        {
            getline(file_below, cur_line);
            if(!cur_line.empty())
            {
                length = cur_line.length();
                VectorXd v(length);
                for(int i = 0; i < length; i++)
                {
                    v[i] = cur_line[i];
                }
                basis.push_back(v);
            }
        }

        file_below.close();

        //fix bug that adds 48 to each entry:
        for(int i = 0; i < basis.size(); i++)
        {
            for(int j = 0; j < length; j++)
            {
                if((basis.at(i)[j] != 1) && (basis.at(i)[j] != 0))
                {
                    basis.at(i)[j] -= 48;
                }
            }
        }

        ofstream new_file;
        new_file.open(to_string(dim_below * 2) + "bit.txt");
        string to_write;

        //goes through each line of the new file
        for(int i = 0; i < (dim_below * 2); i++)
        {
            to_write = "";
            if(i < dim_below)
            {
                for(int j = 0; j < dim_below; j++)
                {
                    if(basis.at(i)[j] == 0)
                    {
                        to_write += "0";
                    }
                    else
                    {
                        to_write += "1";
                    }
                }

                to_write += to_write;
            }

            //last vector
            else if(i == (dim_below * 2) - 1)
            {
                for(int j = 0; j < ((dim_below * 2) - 1); j++)
                {
                    to_write += "0";
                }
                to_write += "1";
            }

            else
            {
                for(int j = 0; j < dim_below; j++)
                {
                    to_write += "0";
                } 

                for(int j = 0; j < dim_below; j++)
                {
                    if(basis.at(i - dim_below)[j] == 0)
                    {
                        to_write += "0";
                    }
                    else
                    {
                        to_write += "1";
                    }
                }  
            }

            //now write this line to the file
            new_file << to_write << "\n";
        }

        new_file.close();
        basis.clear();
        //increment the dimension to the next power of 2
        dim_below = dim_below * 2;
    }
}

// prints the echelon form of some user-entered matrix to the terminal (mostly for testing)
void viewRowReducedMatrix()
{
    vector<VectorXd> basis = getBasis();
    int basis_size = basis.size();
    
    MatrixXd mat(basis_size, basis_size);
    for(int i = 0; i < basis_size; i++)
    {
        mat.col(i) = basis.at(i);
    }

    mat = rowReduce(mat);

    cout << "Here is your row reduced basis matrix: " << endl << endl << mat << endl;

    return;
}

// prints all vectors divided into their subspaces of a user entered dimension
void viewSubspaces()
{
    cout << "What dimension would you like to see the subspaces for?" << endl;
    int dim;
    cin >> dim;
    system("clear");
    while(!((dim != 0) && ((dim & (dim - 1)) == 0)))
    {
        cout << "The dimension you entered is not a power of 2, please try again" << endl << endl << "What dimension would you like to see the subspace for?" << endl;
        cin >> dim;
        system("clear");
    }

    //go through every sequence and compare it to the basis
    vector<vector<VectorXd>> subspaces;

    //get the basis
    ifstream file;
    int length;
    file.open(to_string(dim) + "bit.txt");
    if(file.fail())
    {
        createBasis(dim);
        file.open(dim + "bit.txt");
    }

    string cur_line;
    while(!file.eof())
    {
        getline(file, cur_line);
        if(!cur_line.empty())
        {
            length = cur_line.length();
            vector<VectorXd> v;
            VectorXd newv(length);
            for(int i = 0; i < length; i++)
            {
                newv[i] = cur_line[i];
            }
            v.push_back(newv);
            subspaces.push_back(v);
        }
    }

    file.close();

    //fix bug that adds 48 to each entry:
    for(int i = 0; i < subspaces.size(); i++)
    {
        for(int j = 0; j < length; j++)
        {
            if((subspaces.at(i).at(0)[j] != 1) && (subspaces.at(i).at(0)[j] != 0))
            {
                subspaces.at(i).at(0)[j] -= 48;
            }
        }
    }

    //now add the rest of the sequences
    bool found;
    int sub;
    VectorXd v(dim);
    string sequence;
    for(int i = 1; i <= pow(2, dim); i++)
    {
        printLoadingBar(i, pow(2, dim));
        sequence = decimalToBinary(i);
        while(sequence.length() < dim)
        {
            sequence = "0" + sequence;
        }

        for(int j = 0; j < dim; j++)
        {
            v[j] = (int)sequence[j];
            if((v[j] != 0) && (v[j] != 1))
            {
                v[j] -= 48;
            }
        }

        // if it's odd it goes in the last subspace (this might cut down run-time)
        if(i % 2 != 0)
        {
            if(v != subspaces.at(dim - 1).at(0))
            {
                subspaces.at(dim - 1).push_back(v);
            }
        }

        else
        {
            found = false;
            sub = 0;
            while(!found)
            {
                // this sequence is in the basis
                if(subspaces.at(sub).at(0) == v)
                {
                    break;
                }

                // create a matrix with the basis up to this subspace and our new vector
                MatrixXd mat(dim, sub + 2);
                for(int j = 0; j <= sub; j++)
                {
                    mat.col(j) = subspaces.at(j).at(0);
                }

                mat.col(sub + 1) = v;

                int rank = getRank(mat);
                if(rank == sub + 1)
                {
                    subspaces.at(sub).push_back(v);
                    break;
                }

                else
                {
                    sub++;
                }

            }
        }
    }

    system("clear");

    for(int i = 0; i < subspaces.size(); i++)
    {
        cout << "Here are the vectors for V" << i + 1 << endl;
        // turn them into strings so they print horizontally
        for(int j = 0; j < subspaces.at(i).size(); j++)
        {
            for(int k = 0; k < dim; k++)
            {
                cout << subspaces.at(i).at(j)[k];
            }
            cout << endl;
        }

        cout << endl << endl << endl << endl;
    }

    return;
}




// -----------------------  FUNCTIONS FOR SOLVING CODE JAM ON FLIP MAP   -------------------------

// checks if a subspace with a given matrix is closed under flipping the first two quarters
bool isClosedUnderFlip(vector<VectorXd> basis)
{
    if(!isValidBasis(basis))
    {
        return false;
    }

    int size = basis.size();

    //store rank of matrix before flipping
    MatrixXd mat(basis.at(0).rows(), size);
    for(int i = 0; i < size; i++)
    {
        mat.col(i) = basis.at(i);
    }
    int pre_rank = getRank(mat);

    // make flipped vectors
    // for(int i = 0; i < size; i++)
    // {
    //     VectorXd v(basis.at(0).rows());
    //     VectorXd firsth(basis.at(0).rows() / 4);
    //     VectorXd secondh(basis.at(0).rows() / 4);

    //     for(int j = 0; j < firsth.rows(); j++)
    //     {
    //         firsth[j] = basis.at(i)[j];
    //         secondh[j] = basis.at(i)[j + (basis.at(0).rows() / 4)];
    //     }

    //     for(int j = 0; j < v.rows() / 4; j++)
    //     {
    //         v[j] = secondh[j];
    //         v[j + basis.at(0).rows() / 4] = firsth[j];
    //     }

    //     for(int j = firsth.rows() - 1; j < basis.at(0).rows(); j++)
    //     {
    //         v[j] = basis.at(i)[j];
    //     }

    //     basis.push_back(v);
    // }

    // this code assumes vectors are 4 bits
    for(int i = 0; i < size; i++)
    {
        VectorXd v(4);
        v = basis.at(i);
        int temp = v[0];
        v[0] = v[1];
        v[1] = temp;
        basis.push_back(v);
    }

    MatrixXd post_mat(basis.at(0).rows(), basis.size());
    for(int i = 0; i < basis.size(); i++)
    {
        post_mat.col(i) = basis.at(i);
    }

    basis.clear();

    // check if rank of matrix stayed the same - meaning it's closed under flip
    if(pre_rank == getRank(post_mat))
    {
        return true;
    }

    return false;
}

// checks if 2 2D spaces are the same (not checking chains, just the 2D spaces)
bool areDiff2DSpaces(vector<VectorXd> v1, vector<VectorXd> v2)
{
    // put them in a matrix and check that the rank is not 2
    MatrixXd mat(v1.at(0).rows(), 4);
    mat.col(0) = v1.at(0);
    mat.col(1) = v1.at(1);
    mat.col(2) = v2.at(0);
    mat.col(3) = v2.at(1);

    if(getRank(mat) == 2)
    {
        return false;
    }

    return true;
}

// checks if 2 3D spaces are the same (not checking chains, just the 3D spaces)
bool areDiff3DSpaces(vector<VectorXd> v1, vector<VectorXd> v2)
{
    MatrixXd mat(v1.at(0).rows(), 6);
    mat.col(0) = v1.at(0);
    mat.col(1) = v1.at(1);
    mat.col(2) = v1.at(2);
    mat.col(3) = v2.at(0);
    mat.col(4) = v2.at(1);
    mat.col(5) = v2.at(2);

    if(getRank(mat) == 3)
    {
        return false;
    }

    return true;
}

// finds all the different 2D invariant subspaces and prints them to a txt file
void twosAreDiffUnderFlip()
{
    // make a vector of all the sets of 2 vectors that are invariant, then go through them all and check that they are different
    
    vector<vector<VectorXd>> found_spaces;

    vector<VectorXd> first_dim;
    VectorXd v(4);
    v << 1, 1, 1, 1;
    first_dim.push_back(v);
    v << 1, 1, 1, 0;
    first_dim.push_back(v);
    v << 1, 1, 0, 1;
    first_dim.push_back(v);
    v << 1, 1, 0, 0;
    first_dim.push_back(v);
    v << 0, 0, 1, 1;
    first_dim.push_back(v);
    v << 0, 0, 1, 0;
    first_dim.push_back(v);
    v << 0, 0, 0, 1;
    first_dim.push_back(v);

    for(int i = 0; i < 7; i++)
    {
        vector<VectorXd> cur_space;
        cur_space.push_back(first_dim.at(i));
        for(int j = 1; j < 16; j++)
        {
            string bin = decimalToBinary(j);
            while(bin.length() < 4)
            {
                bin = "0" + bin;
            }

            VectorXd second(4);
            for(int k = 0; k < 4; k++)
            {
                second[k] = bin[k];
                if((second[k] != 1) && (second[k] != 0))
                {
                    second[k] -= 48;
                }
            }

            cur_space.push_back(second);
            // check if this is an invariant space, then check that it's not the same as spaces already in found_spaces
            if(!isClosedUnderFlip(cur_space))
            {
                cur_space.pop_back();
                continue;
            }


            bool exists = false;

            for(int p = 0; p < found_spaces.size(); p++)
            {
                // cout << "First Basis: " << found_spaces.at(p).at(0)[0] << found_spaces.at(p).at(0)[1] << found_spaces.at(p).at(0)[2] << found_spaces.at(p).at(0)[3] << "    " << found_spaces.at(p).at(1)[0] << found_spaces.at(p).at(1)[1] << found_spaces.at(p).at(1)[2] << found_spaces.at(p).at(1)[3] << endl;
                // cout << "Second Basis: " << cur_space.at(0)[0] << cur_space.at(0)[1] << cur_space.at(0)[2] << cur_space.at(0)[3] << "    " << cur_space.at(1)[0] << cur_space.at(1)[1] << cur_space.at(1)[2] << cur_space.at(1)[3] << endl;
                if(!areDiff2DSpaces(found_spaces.at(p), cur_space))
                {
                    exists = true;
                    break;
                }
            }

            // 1100 0011 and 0011 1100 are both showing up as different spaces and not going into the if statement
            // the function is mistakenly saying that they are different spaces when they are obviously not. fix this

            if(!exists)
            {
                found_spaces.push_back(cur_space);
            }

            cur_space.pop_back();
        }

        cur_space.pop_back();
    }

    ofstream file;
    file.open("twos_under_flip.txt");
    file << "Here are the different 2D spaces under flip:" << endl << endl;
    for(int i = 0; i < found_spaces.size(); i++)
    {
        file << found_spaces.at(i).at(0)[0] << found_spaces.at(i).at(0)[1] << found_spaces.at(i).at(0)[2] << found_spaces.at(i).at(0)[3] << "      " << found_spaces.at(i).at(1)[0] << found_spaces.at(i).at(1)[1] << found_spaces.at(i).at(1)[2] << found_spaces.at(i).at(1)[3] << endl;
    }

    file.close();


    return;
}

// finds and writes to "flip_matches.txt" all of the 1D spaces matched with which 2D spaces they fit in
void findWhichTwosUnderFlip()
{
    ifstream reading;
    reading.open("twos_under_flip.txt");
    string line;

    // so it skips the first 2 lines
    getline(reading, line);
    getline(reading, line);
    vector<vector<VectorXd>> twos;

    // get all the unique 2D spaces
    for(int i = 0; i < 11; i++)
    {
        vector<VectorXd> newvs;
        VectorXd v(4);
        getline(reading, line);

        v << line[0], line[1], line[2], line[3];
        for(int j = 0; j < 4; j++)
        {
            if((v[j] != 0) && (v[j] != 1))
            {
                v[j] -= 48;
            }
        }
        newvs.push_back(v);

        v << line[10], line[11], line[12], line[13];
        for(int j = 0; j < 4; j++)
        {
            if((v[j] != 0) && (v[j] != 1))
            {
                v[j] -= 48;
            }
        }
        newvs.push_back(v);

        twos.push_back(newvs);
    }

    reading.close();

    // find which 2D spaces the 1D spaces match to, put them in a file
    vector<VectorXd> ones;
    VectorXd v(4);

    v << 1, 1, 1, 1;
    ones.push_back(v);
    v << 1, 1, 0, 0;
    ones.push_back(v);
    v << 1, 1, 1, 0;
    ones.push_back(v);
    v << 1, 1, 0, 1;
    ones.push_back(v);
    v << 0, 0, 1, 1;
    ones.push_back(v);
    v << 0, 0, 0, 1;
    ones.push_back(v);
    v << 0, 0, 1, 0;
    ones.push_back(v);

    ofstream writing;
    writing.open("flip_matches.txt");
    writing << "These are all the matches of 1D spaces to unique 2D spaces:" << endl << endl;

    for(int i = 0; i < 7; i++)
    {
        MatrixXd mat(4, 3);
        mat.col(0) = ones.at(i);
        for(int j = 0; j < twos.size(); j++)
        {
            mat.col(1) = twos.at(j).at(0);
            mat.col(2) = twos.at(j).at(1);
            if(getRank(mat) == 2)
            {
                writing << ones.at(i)[0] << ones.at(i)[1] << ones.at(i)[2] << ones.at(i)[3] << "      " << twos.at(j).at(0)[0] << twos.at(j).at(0)[1] << twos.at(j).at(0)[2] << twos.at(j).at(0)[3] << "    " << twos.at(j).at(1)[0] << twos.at(j).at(1)[1] << twos.at(j).at(1)[2] << twos.at(j).at(1)[3] << endl;
            }
        }

        writing << endl;
    }

    writing.close();

    return;
}

// write this ^ backwards
void findWhichOnesUnderFlip()
{
    ifstream reading;
    reading.open("twos_under_flip.txt");
    string line;

    // so it skips the first 2 lines
    getline(reading, line);
    getline(reading, line);
    vector<vector<VectorXd>> twos;

    // get all the unique 2D spaces
    for(int i = 0; i < 11; i++)
    {
        vector<VectorXd> newvs;
        VectorXd v(4);
        getline(reading, line);

        v << line[0], line[1], line[2], line[3];
        for(int j = 0; j < 4; j++)
        {
            if((v[j] != 0) && (v[j] != 1))
            {
                v[j] -= 48;
            }
        }
        newvs.push_back(v);

        v << line[10], line[11], line[12], line[13];
        for(int j = 0; j < 4; j++)
        {
            if((v[j] != 0) && (v[j] != 1))
            {
                v[j] -= 48;
            }
        }
        newvs.push_back(v);

        twos.push_back(newvs);
    }

    reading.close();

    // find which 1D spaces the 2D spaces match to, put them in a file
    vector<VectorXd> ones;
    VectorXd v(4);

    v << 1, 1, 1, 1;
    ones.push_back(v);
    v << 1, 1, 0, 0;
    ones.push_back(v);
    v << 1, 1, 1, 0;
    ones.push_back(v);
    v << 1, 1, 0, 1;
    ones.push_back(v);
    v << 0, 0, 1, 1;
    ones.push_back(v);
    v << 0, 0, 0, 1;
    ones.push_back(v);
    v << 0, 0, 1, 0;
    ones.push_back(v);

    ofstream writing;
    writing.open("flip_matches_backwards.txt");
    writing << "These are all the matches of 2D spaces to unique 1D spaces:" << endl << endl;

    for(int i = 0; i < twos.size(); i++)
    {
        MatrixXd mat(4, 3);
        mat.col(0) = twos.at(i).at(0);
        mat.col(1) = twos.at(i).at(1);
        for(int j = 0; j < 7; j++)
        {
            mat.col(2) = ones.at(j);
            if(getRank(mat) == 2)
            {
                writing << twos.at(i).at(0)[0] << twos.at(i).at(0)[1] << twos.at(i).at(0)[2] << twos.at(i).at(0)[3] << "    " << twos.at(i).at(1)[0] << twos.at(i).at(1)[1] << twos.at(i).at(1)[2] << twos.at(i).at(1)[3] << "        " << ones.at(j)[0] << ones.at(j)[1] << ones.at(j)[2] << ones.at(j)[3] << endl;
            }
        }

        writing << endl;
    }

    writing.close();

    return;
}

// find all the unique 3D spaces
void findDiffThreesUnderFlip()
{
    // make a vector of all the sets of 2 vectors that are invariant, then go through them all and check that they are different
    
    vector<vector<VectorXd>> found_spaces;

    vector<VectorXd> first_dim;
    VectorXd v(4);
    v << 1, 1, 1, 1;
    first_dim.push_back(v);
    v << 1, 1, 1, 0;
    first_dim.push_back(v);
    v << 1, 1, 0, 1;
    first_dim.push_back(v);
    v << 1, 1, 0, 0;
    first_dim.push_back(v);
    v << 0, 0, 1, 1;
    first_dim.push_back(v);
    v << 0, 0, 1, 0;
    first_dim.push_back(v);
    v << 0, 0, 0, 1;
    first_dim.push_back(v);

    for(int i = 0; i < 7; i++)
    {
        vector<VectorXd> cur_space;
        cur_space.push_back(first_dim.at(i));
        for(int j = 1; j < 16; j++)
        {
            string bin = decimalToBinary(j);
            while(bin.length() < 4)
            {
                bin = "0" + bin;
            }

            VectorXd second(4);
            for(int k = 0; k < 4; k++)
            {
                second[k] = bin[k];
                if((second[k] != 1) && (second[k] != 0))
                {
                    second[k] -= 48;
                }
            }

            cur_space.push_back(second);
            // check if this is an invariant space, then check that it's not the same as spaces already in found_spaces
            if(!isClosedUnderFlip(cur_space))
            {
                cur_space.pop_back();
                continue;
            }

            // find 3rd vec
            for(int k = 1; k < 16; k++)
            {
                string binary = decimalToBinary(k);
                while(binary.length() < 4)
                {
                    binary = "0" + binary;
                }

                VectorXd third(4);
                for(int h = 0; h < 4; h++)
                {
                    third[h] = binary[h];
                    if((third[h] != 0) && (third[h] != 1))
                    {
                        third[h] -= 48;
                    }
                }

                cur_space.push_back(third);

                // check it's still invariant
                if(!isClosedUnderFlip(cur_space))
                {
                    cur_space.pop_back();
                    continue;
                }

                // check it's not already in the found_spaces
                bool exists = false;
                for(int p = 0; p < found_spaces.size(); p++)
                {
                    if(!areDiff3DSpaces(found_spaces.at(p), cur_space))
                    {
                        exists = true;
                        break;
                    }
                }

                if(!exists)
                {
                    found_spaces.push_back(cur_space);
                }

                cur_space.pop_back();
            }

            cur_space.pop_back();
        }

        cur_space.pop_back();
    }

    ofstream file;
    file.open("threes_under_flip.txt");
    file << "Here are the unique 3D spaces under flip:" << endl << endl;
    for(int i = 0; i < found_spaces.size(); i++)
    {
        file << found_spaces.at(i).at(0)[0] << found_spaces.at(i).at(0)[1] << found_spaces.at(i).at(0)[2] << found_spaces.at(i).at(0)[3] << "      " << found_spaces.at(i).at(1)[0] << found_spaces.at(i).at(1)[1] << found_spaces.at(i).at(1)[2] << found_spaces.at(i).at(1)[3] << "      " << found_spaces.at(i).at(2)[0] << found_spaces.at(i).at(2)[1] << found_spaces.at(i).at(2)[2] << found_spaces.at(i).at(2)[3] << endl;;
    }

    file.close();


    return;
}




// -----------------------  FUNCTIONS FOR SOLVING CODE JAM ON FLIP MAP   -------------------------

// checks if a subspace with a given matrix is closed under swapping the halves of the vectors
bool isClosedUnderSwap(vector<VectorXd> basis)
{
    if(!isValidBasis(basis))
    {
        return false;
    }

    int size = basis.size();

    //store rank of matrix before flipping
    MatrixXd mat(basis.at(0).rows(), size);
    for(int i = 0; i < size; i++)
    {
        mat.col(i) = basis.at(i);
    }

    int pre_rank = getRank(mat);

    // make swapped vectors
    for(int i = 0; i < size; i++)
    {
        VectorXd first(basis.at(0).rows() / 2);
        for(int j = 0; j < first.rows(); j++)
        {
            first[j] = basis.at(i)[j];
        }

        VectorXd v(basis.at(0).rows());
        for(int j = 0; j < v.rows() / 2; j++)
        {
            v[j] = basis.at(i)[j + (v.rows() / 2)];
            v[j + (v.rows() / 2)] = first[j];
        }

        basis.push_back(v);
    }

    MatrixXd post_mat(basis.at(0).rows(), basis.size());
    for(int i = 0; i < basis.size(); i++)
    {
        post_mat.col(i) = basis.at(i);
    }

    basis.clear();

    // check if rank of matrix stayed the same - meaning it's closed under swap
    if(pre_rank == getRank(post_mat))
    {
        return true;
    }

    return false;
}

// finds all chains for the map that swaps the first two bits and prints them to a txt file
void fourBitChains()
{
    // holds all of the bases we find
    vector<vector<VectorXd>> found_spaces;

    // holds all of the vectors that are a basis for a 1d invariant space
    vector<VectorXd> temp;
    VectorXd v(4);
    v << 1, 1, 0, 0;
    temp.push_back(v);

    v << 1, 1, 1, 1;
    temp.push_back(v);

    v << 1, 1, 1, 0;
    temp.push_back(v);

    v << 1, 1, 0, 1;
    temp.push_back(v);

    v << 0, 0, 1, 1;
    temp.push_back(v);

    v << 0, 0, 1, 0;
    temp.push_back(v);

    v << 0, 0, 0, 1;
    temp.push_back(v);

    // loop through each of these first vectors
    for(int p = 0; p < 7; p++)
    {
        // the_vec stores the temporary basis that we are working on creating, once it is a complete basis, it is added to found_spaces and tries to create another
        vector<VectorXd> the_vec;
        the_vec.push_back(temp.at(p));

        for(int i = 1; i < 16; i++)
        {
            // create the string of the next vec
            string bin = decimalToBinary(i);
            while(bin.length() < 4)
            {
                bin = "0" + bin;
            }

            // turn the string into a vector
            VectorXd v1(4);
            for(int j = 0; j < 4; j++)
            {
                v1[j] = bin[j];
                if((v1[j] != 0) && (v1[j] != 1))
                {
                    v1[j] -= 48;
                }
            }

            the_vec.push_back(v1);

            // if this new basis is not an invariant basis by our definition, scrap it and go to the next vector
            if(!isClosedUnderFlip(the_vec))
            {
                the_vec.pop_back();
                continue;
            }

            // start trying to add another vec
            for(int j = 1; j < 16; j++)
            {
                // create new vec in string form
                bin = decimalToBinary(j);
                while(bin.length() < 4)
                {
                    bin = "0" + bin;
                }

                // turn string into actual vector
                for(int k = 0; k < 4; k++)
                {
                    v1[k] = bin[k];
                    if((v1[k] != 0) && (v1[k] != 1))
                    {
                        v1[k] -= 48;
                    }
                }

                the_vec.push_back(v1);

                // if this creates a non-invariant subspace, scrap it and check the next vector
                if(!isClosedUnderFlip(the_vec))
                {
                    the_vec.pop_back();
                    continue;
                }

                // start trying to add the next vector
                for(int k = 1; k < 15; k++)
                {
                    // make a vec in string form
                    bin = decimalToBinary(k);
                    while(bin.length() < 4)
                    {
                        bin = "0" + bin;
                    }

                    // turn the string into an actual vector
                    for(int h = 0; h < 4; h++)
                    {
                        v1[h] = bin[h];
                        if((v1[h] != 0) && (v1[h] != 1))
                        {
                            v1[h] -= 48;
                        }
                    }

                    the_vec.push_back(v1);

                    // if this creates a non-invariant subspace, scrap this vector and move to the next one
                    if(!isClosedUnderFlip(the_vec))
                    {
                        the_vec.pop_back();
                        continue;
                    }

                    found_spaces.push_back(the_vec);
                    the_vec.pop_back();
                }

                the_vec.pop_back();
            }

            the_vec.pop_back();
        }
        the_vec.pop_back();
    }
    
    // this makes a vector of bases for the subspaces, however they have not been worked through to eliminate bases that cover the same spaces
    // to check this, put the bases in a matrix and check that it only has a rank of 4

    cout << "size = " << found_spaces.size() << endl;
    for(int i = 0; i < found_spaces.size(); i++)
    {
        // cout << "size = " << found_spaces.size() << endl;
        for(int j = i + 1; j < found_spaces.size(); j++)
        {
            if(!areDiffSpaces(found_spaces.at(i), found_spaces.at(j)))
            {
                // then all the subspaces are the same so we should erase one of the bases
                found_spaces.erase(found_spaces.begin() + j);
                j--;
            }
        }
    }

    // now everything in vecs should be bases for different chains
    // let's print them
    ofstream file;
    file.open("flip.txt");
    file << "There are " << found_spaces.size() << " different chains of subspaces for sigma where sigma switches the first two quarters of the vector." << endl << "Here they are:" << endl << endl;
    for(int i = 0; i < found_spaces.size(); i++)
    {
        for(int j = 0; j < 4; j++)
        {
            file << found_spaces.at(i).at(j)[0] << found_spaces.at(i).at(j)[1] << found_spaces.at(i).at(j)[2] << found_spaces.at(i).at(j)[3] << "     ";
        }
        file << endl;
    }

    found_spaces.clear();
    temp.clear();

    file.close();
}

// finds all chains for the sigma that swaps the halves and prints them to a txt file
void fourBitHalfSwapChains()
{
    // holds all of the bases we find
    vector<vector<VectorXd>> found_spaces;

    // holds all of the vectors that are a basis for a 1d invariant space
    vector<VectorXd> temp;
    VectorXd v(4);
    v << 1, 1, 1, 1;
    temp.push_back(v);

    v << 1, 0, 1, 0;
    temp.push_back(v);

    v << 0, 1, 0, 1;
    temp.push_back(v);

    // loop through each of these first vectors
    for(int p = 0; p < 3; p++)
    {
        // the_vec stores the temporary basis that we are working on creating, once it is a complete basis, it is added to found_spaces and tries to create another
        vector<VectorXd> the_vec;
        the_vec.push_back(temp.at(p));

        for(int i = 1; i < 16; i++)
        {
            // create the string of the next vec
            string bin = decimalToBinary(i);
            while(bin.length() < 4)
            {
                bin = "0" + bin;
            }

            // turn the string into a vector
            VectorXd v1(4);
            for(int j = 0; j < 4; j++)
            {
                v1[j] = bin[j];
                if((v1[j] != 0) && (v1[j] != 1))
                {
                    v1[j] -= 48;
                }
            }

            the_vec.push_back(v1);

            // if this new basis is not an invariant basis by our definition, scrap it and go to the next vector
            if(!isClosedUnderSwap(the_vec))
            {
                the_vec.pop_back();
                continue;
            }

            // start trying to add another vec
            for(int j = 1; j < 16; j++)
            {
                // create new vec in string form
                bin = decimalToBinary(j);
                while(bin.length() < 4)
                {
                    bin = "0" + bin;
                }

                // turn string into actual vector
                for(int k = 0; k < 4; k++)
                {
                    v1[k] = bin[k];
                    if((v1[k] != 0) && (v1[k] != 1))
                    {
                        v1[k] -= 48;
                    }
                }

                the_vec.push_back(v1);

                // if this creates a non-invariant subspace, scrap it and check the next vector
                if(!isClosedUnderSwap(the_vec))
                {
                    the_vec.pop_back();
                    continue;
                }

                // start trying to add the next vector
                for(int k = 1; k < 15; k++)
                {
                    // make a vec in string form
                    bin = decimalToBinary(k);
                    while(bin.length() < 4)
                    {
                        bin = "0" + bin;
                    }

                    // turn the string into an actual vector
                    for(int h = 0; h < 4; h++)
                    {
                        v1[h] = bin[h];
                        if((v1[h] != 0) && (v1[h] != 1))
                        {
                            v1[h] -= 48;
                        }
                    }

                    the_vec.push_back(v1);

                    // if this creates a non-invariant subspace, scrap this vector and move to the next one
                    if(!isClosedUnderSwap(the_vec))
                    {
                        the_vec.pop_back();
                        continue;
                    }

                    found_spaces.push_back(the_vec);
                    the_vec.pop_back();
                }

                the_vec.pop_back();
            }

            the_vec.pop_back();
        }
        the_vec.pop_back();
    }
    
    // this makes a vector of bases for the subspaces, however they have not been worked through to eliminate bases that cover the same spaces
    // to check this, put the bases in a matrix and check that it only has a rank of 4

    cout << "size = " << found_spaces.size() << endl;
    for(int i = 0; i < found_spaces.size(); i++)
    {
        // cout << "size = " << found_spaces.size() << endl;
        for(int j = i + 1; j < found_spaces.size(); j++)
        {
           
            if(!areDiffSpaces(found_spaces.at(i), found_spaces.at(j)))
            {
                // then all the subspaces are the same so we should erase one of the bases
                found_spaces.erase(found_spaces.begin() + j);
                j--;
            }
        }
    }

    // now everything in vecs should be bases for different chains
    // let's print them
    ofstream file;
    file.open("swap.txt");
    file << "There are " << found_spaces.size() << " different chains of subspaces for sigma where sigma switches the halves of the vectors." << endl << "Here they are:" << endl << endl;
    for(int i = 0; i < found_spaces.size(); i++)
    {
        for(int j = 0; j < 4; j++)
        {
            file << found_spaces.at(i).at(j)[0] << found_spaces.at(i).at(j)[1] << found_spaces.at(i).at(j)[2] << found_spaces.at(i).at(j)[3] << "     ";
        }
        file << endl;
    }

    found_spaces.clear();
    temp.clear();

    file.close();
}




// -------------------------------------  TESTS AND MAIN   ---------------------------------------

// test functions here
void tests()
{
    findDiffThreesUnderFlip();
    findWhichTwosUnderFlip();
    twosAreDiffUnderFlip();
    findWhichOnesUnderFlip();

    vector<VectorXd> v1;
    vector<VectorXd> v2;
    VectorXd v(4);
    v << 1, 1, 0, 0;
    v1.push_back(v);
    v << 0, 0, 1, 1;
    v2.push_back(v);
    v2.push_back(v1.at(0));
    v1.push_back(v);

    // cout << "the bases are different: " << areDiffSpaces(v1, v2) << endl;

    // vector<VectorXd> v1;
    // VectorXd v(4);
    // v << 1,1,1,1;
    // v1.push_back(v);
    // v << 1,1,0,0;
    // v1.push_back(v);
    // v << 0,0,0,1;
    // v1.push_back(v);
    // v << 1,0,1,1;
    // v1.push_back(v);

    // isClosedUnderSwap(v1);

    // vector<VectorXd> v2 = v1;


    // bool work = areDiffSpaces(v1, v2);
    // if(work)
    // {
    //     cout << "Uh oh, function doesn't work." << endl;
    // }

    // else
    // {
    //     cout << "Yay it works!" << endl;
    // }
}

int main()
{
    system("clear");
    tests();
    string input = "0";
    while(input != "7")
    {
        cout << "Choose what you would like to do:" << endl;
        cout << "1. Create a new basis file" << endl << "2. Check the validity of an existing basis" << endl << "3. View the row reduced form of a basis matrix"  << endl << "4. View the subspaces of a dimension of sequences" << endl << "5. View the chains for 4 bits for flip sigma" << endl << "6. View the chains for 4 bits for swap sigma" << endl << "7. Exit the program" << endl;
        cin >> input;
        system("clear");
        while ((input != "1") && (input != "2") && (input != "3") && (input != "4") && (input != "5") && (input != "6") && (input != "7"))
        {
            cout << "Invalid input. Please try again and pick one of the listed options:" << endl;
            cout << "1. Create a new basis file" << endl << "2. Check the validity of an existing basis" << endl << "3. View the row reduced form of a basis matrix" << endl << "4. View the subspaces of a dimension of sequences" << endl << "5. View the chains for 4 bits for special sigma" << endl << "6. View the chains for 4 bits for swap sigma" << endl << "7. Exit the program" << endl;
            cin >> input;
            system("clear");
        }

        if(input == "1")
        {
            cout << "What dimension of basis would you like?" << endl;
            int dim;
            cin >> dim;
            system("clear");
            createBasis(dim);
        }

        else if(input == "2")
        {
            vector<VectorXd> basis = getBasis();
            bool valid = isValidBasis(basis);
            bool invariant = isInvariantBasis(basis);
            basis.clear();

            if(valid && invariant)
            {
                cout << "The vectors you entered form an invariant basis." << endl;
            }

            else if(valid)
            {
                cout << "The vectors you entered form a non-invariant basis." << endl;
            }

            else
            {
                cout << "The vectors you entered do not form a basis." << endl; 
            }
        }

        else if(input == "3")
        {
            viewRowReducedMatrix();
        }

        else if(input == "4")
        {
            viewSubspaces();
        }

        else if(input == "5")
        {
            fourBitChains();
        }

        else if(input == "6")
        {
            fourBitHalfSwapChains();
        }

        else
        {
            cout << "Exiting program now." << endl;
            return 0;
        }

        cout << endl;
    }
    
    return 0;
}