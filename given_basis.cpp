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
    bool pivot_exists;
    int starting_row;

    /** Algorithm:
     * find a row that contains a one in the current column
     * swap ith row, where i is the current column
     * add this row to every row below it that contains a one in the ith column
    */

   for(int i = 0; i < cols; i++)
   {
        pivot_exists = false;

        //this loop accounts for if the previous column was not a pivot, then the function needs to place the next pivot in the row of the previous
        starting_row = i;
        while((i != 0) && (mat(starting_row - 1, i-1) != 1))
        {
            starting_row--;
        }
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

    // clear the vector for memory purposes before returning
    basis.clear();

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

    //we already know we have a basis, for it to be invariant, the rank of the matrix must be the same with the cycled vectors as without
    if(rank == basis_size)
    {
        return true;
    }

    return false;
}

//function that will create a basis for an inputted dimension and put it in a file
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
     * see if we have a file for the dimension down - if not, we're gonna make one
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
            first_file << "1010" << "\n";
            first_file << "1001" << "\n";
            first_file << "1000" << "\n";
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
                to_write = "1";
                for(int j = 0; j < ((dim_below * 2) - 1); j++)
                {
                    to_write += "0";
                }
            }

            else
            {
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

                for(int j = 0; j < dim_below; j++)
                {
                    to_write += "0";
                }   
            }

            //now write this line to the file
            new_file << to_write << "\n";
        }

        new_file.close();
        //increment the dimension to the next power of 2
        dim_below = dim_below * 2;
    }
}

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

int main()
{
    string input = "0";
    while(input != "4")
    {
        cout << "Choose what you would like to do:" << endl;
        cout << "1. Create a new basis file" << endl << "2. Check the validity of an existing basis" << endl << "3. View the row reduced form of a basis matrix"  << endl << "4. Exit the program" << endl;
        cin >> input;
        while ((input != "1") && (input != "2") && (input != "3") && (input != "4"))
        {
            cout << "Invalid input. Please try again and pick one of the listed options:" << endl;
            cout << "1. Create a new basis file" << endl << "2. Check the validity of an existing basis" << endl << "3. View the row reduced form of a basis matrix"  << endl << "4. Exit the program" << endl;
            cin >> input;
        }

        if(input == "1")
        {
            cout << "What dimension of basis would you like?" << endl;
            int dim;
            cin >> dim;
            createBasis(dim);
        }

        else if(input == "2")
        {
            vector<VectorXd> basis = getBasis();
            bool valid = isValidBasis(basis);
            bool invariant = isInvariantBasis(basis);

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

        else
        {
            cout << "Exiting program now." << endl;
            return 0;
        }

        cout << endl;
    }
    
    return 0;
}