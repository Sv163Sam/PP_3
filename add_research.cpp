#include "add.hpp"

void _create_file()
{
    int row = 0, column = 0;
    int** matrix;
    
    srand((unsigned)time(NULL));
    
    row = 10;
    column = 10;
    
    matrix = new int* [row];
    for(size_t i = 0; i < row; i++)
    {
        matrix[i] = new int [column];
        for(size_t j = 0; j < column; j++)
        {
            matrix[i][j] = rand() % 10;
        }
    }
    
    FILE* out = fopen("./input_f.txt", "wt");
    fprintf(out, "%d", row);
    fprintf(out, " ");
    fprintf(out, "%d", column);
    for (size_t i = 0; i < row; i++)
    {
        fprintf(out, "\n");
        for(size_t j = 0; j < column; j++)
        {
            fprintf(out, "%d", matrix[i][j]);
            fprintf(out, " ");
        }
    }
    fclose(out);
    
    int tmp = row;
    row = column;
    column = tmp;
    
    matrix = new int* [row];
    for(size_t i = 0; i < row; i++)
    {
        matrix[i] = new int [column];
        for(size_t j = 0; j < column; j++)
        {
            matrix[i][j] = rand() % 10;
        }
    }
    
    out = fopen("./input_s.txt", "wt");
    fprintf(out, "%d", row);
    fprintf(out, " ");
    fprintf(out, "%d", column);
    for (size_t i = 0; i < row; i++)
    {
        fprintf(out, "\n");
        for(size_t j = 0; j < column; j++)
        {
            fprintf(out, "%d", matrix[i][j]);
            fprintf(out, " ");
        }
    }
    fclose(out);
}

