#include <ctime>
#include "add.hpp"
#include "mpi.h"

struct mtrx
{
    int** matrix;
    int row_size = 0;
    int column_size = 0;
};

void _print_matrix_(mtrx obj)
{
    std::cout << "Matrix from file" << std::endl;
    std::cout << "Amount of rows - " << obj.row_size << std::endl;
    std::cout << "Amount of columns - " << obj.column_size << std::endl;
    std::cout << "Matrix: " << std::endl;
    for (size_t i = 0; i < obj.row_size; i++)
    {
        std::cout << std::endl;
        for(size_t j = 0; j < obj.column_size; j++)
        {
            std::cout << obj.matrix[i][j] << "\t" << std::endl;
        }
    }
    return;
}

mtrx _initialize_matrix_(const char* path, const char* mode)
{
    mtrx obj;
    int row = 0, column = 0;
    FILE *in = fopen(path, mode);
    
    fscanf(in, "%d", &row);
    fscanf(in, "%d", &column);
    
    int** first_m = new int* [row];
    for(size_t i = 0; i < row; i++)
    {
        first_m[i] = new int [column];
        for(size_t j = 0; j < column; j++)
        {
            fscanf(in, "%d", &first_m[i][j]);
        }
    }
    obj.matrix = first_m;
    obj.row_size = row;
    obj.column_size = column;
    fclose(in);
    return obj;
}

mtrx _multi_matrix_(int *&A, int *&B, int *&C, int &size)
{
    int proc_num = size / 2;
    int size_dupl = size;
    int index;
    double temp = 0;
    
    MPI_Status Status;
    int proc_part_size = size_dupl / proc_num;
    int proc_part_elem = proc_part_size * size_dupl;
    int proc_part = size_dupl / proc_num;
    int part = proc_part * size_dupl;
    
    int* buf_a = new int[proc_part_elem];
    int* buf_b = new int[proc_part_elem];
    int* buf_c = new int[proc_part_elem];
    
    if (proc_rank == 0)
    {
        Flip(B, size);
    }
    
    MPI_Scatter(A, part, MPI_INT, buf_a, part, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, part, MPI_INT, buf_b, part, MPI_INT, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i < proc_part_size; i++)
    {
        for (int j = 0; j < proc_part_size; j++)
        {
            for (int k = 0; k < size_dupl; k++)
            {
                temp += buf_a[i * size_dupl + k] * buf_b[j * size_dupl + k];
            }
            buf_c[i * size_dupl + j + proc_part_size * proc_rank] = temp;
            temp = 0.0;
        }
    }

    int next_proc;
    int prev_proc;
    
    for (int p = 1; p < proc_num; p++)
    {
        next_proc = proc_rank + 1;
        if (proc_rank == proc_num - 1)
            next_proc = 0;
        
        prev_proc = proc_rank - 1;
        
        if (proc_rank == 0)
            prev_proc = proc_num - 1;
        
        MPI_Sendrecv_replace(buf_b, part, MPI_INT, next_proc, 0, prev_proc, 0, MPI_COMM_WORLD, &Status);

        temp = 0.0;
        
        for (int i = 0; i < proc_part_size; i++)
        {
            for (int j = 0; j < proc_part_size; j++)
            {
                for (int k = 0; k < size_dupl; k++)
                {
                    temp += buf_a[i * size_dupl + k] * buf_b[j * size_dupl + k];
                }
                if (proc_rank - p >= 0)
                    index = proc_rank - p;
                
                else
                    index = (proc_num - p + proc_rank);
                
                buf_c[i * size_dupl + j + index * proc_part_size] = temp;
                temp = 0.0;
            }
        }
    }
    
    MPI_Gather(buf_c, proc_part_elem, MPI_INT, C, proc_part_elem, MPI_INT, 0, MPI_COMM_WORLD);

    delete []bufA;
    delete []bufB;
    delete []bufC;
}


/*
mtrx _multi_matrix_(mtrx matrix_f, mtrx matrix_s)
{
    mtrx obj;
    if (matrix_f.column_size != matrix_s.row_size) throw std::length_error("Multiplication is impossible");
    int** res_matrix = new int* [matrix_f.row_size];
    
    int threads;
        omp_set_num_threads(8);
    #pragma omp parallel shared(threads)
    {
        threads = omp_get_num_threads();
#pragma omp for
        for (size_t i = 0; i < matrix_f.row_size; i++)
        {
            res_matrix[i] = new int[matrix_s.column_size];
            for (size_t j = 0; j < matrix_s.column_size; j++)
            {
                res_matrix[i][j] = 0;
                for (size_t k = 0; k < matrix_f.column_size; k++)
                {
                    res_matrix[i][j] += matrix_f.matrix[i][k] * matrix_s.matrix[k][j];
                }
            }
        }
    }
    obj.matrix = res_matrix;
    obj.row_size = matrix_f.row_size;
    obj.column_size = matrix_s.column_size;
    return obj;
}
*/

void _write_file_matrix_(mtrx obj, const char* path, const char* mode, clock_t _time)
{
    FILE* out = fopen(path, mode);
    for (size_t i = 0; i < obj.row_size; i++)
    {
        for(size_t j = 0; j < obj.column_size; j++)
        {
            fprintf(out, "%d", obj.matrix[i][j]);
            fprintf(out, "\t");
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\nRuntime: %lf ms\n", (float)_time/1000);
    fprintf(out, "Amount of rows: %d\n", obj.row_size);
    fprintf(out, "Amount of columns: %d", obj.column_size);
    fclose(out);
    return;
}

int main(int argc, char* argv[])
{
    const char* path_f = "./input_f.txt";
    const char* path_s = "./input_s.txt";
    const char* r_m = "rt";
    
    const char* e_path = "./output.txt";
    const char* w_m = "wt";

    int proc_rank;
    int proc_count;
    
    mtrx first_mtrx;
    mtrx second_mtrx;
    mtrx res_matrix;
    
    _create_file();
    
    first_mtrx = _initialize_matrix_(path_f, r_m);
    second_mtrx = _initialize_matrix_(path_s, r_m);
    
    int** res = new int* [first_mtrx.row_size];
    for(size_t i = 0; i < first_mtrx.row_size; i++)
    {
        res[i] = new int [first_mtrx.column_size];
        for(size_t j = 0; j < first_mtrx.column_size; j++)
        {
            fscanf(in, "%d", &res[i][j]);
        }
    }
    
    res_matrix.matrix = res;
    res_matrix.row_size = first_mtrx.row_size;
    res_matrix.column_size = first_mtrx.column_size;
    
    try
    {
        MPI_INIT(&argc, &argv);
        //MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
        //MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        clock_t start_time = clock();
        _multi_matrix_(first_mtrx.matrix, second_mtrx.matrix, res_matrix.matrix, res_matrix.column_size);//тут сам поставь & * - чо надо будет для передачи аргументов!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //
        //
        //
        //
        //
        clock_t end_time = clock();
        MPI_Finalize();
        
        clock_t time = end_time - start_time;
        std::cout << time << std::endl;
        
        std::cout << "\nC++ calculations: " << std::endl;
        for (size_t i = 0; i < res_matrix.row_size; i++)
        {
            std::cout << std::endl;
            for (size_t j = 0; j < res_matrix.column_size; j++)
            {
                std::cout << res_matrix.matrix[i][j] << "\t";
            }
        }
        _write_file_matrix_(res_matrix, e_path, w_m, time);
    }
    catch (const std::length_error &obj)
    {
        std::cout << "Emergency program crash: " << obj.what() << "\nPlease change the data in the file" << std::endl;
        exit(1);
    }
    std::cout << "\nPython check script\n" << std::endl;

    system("python3 /Users/vladimirskobcov/PycharmProjects/pythonProject5/main.py");
    return 0;
}

