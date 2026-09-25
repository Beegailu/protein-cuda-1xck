#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_ATOMS 60000
#define MAX_RESIDUES 10000
#define CUTOFF 6.0

typedef struct {
    char atom_name[5];
    char res_name[4];
    char chain;
    int res_seq;
    char icode;
    char element[3];

    float x;
    float y;
    float z;

    int residue_index;
} Atom;

typedef struct {
    char res_name[4];
    char chain;
    int res_seq;
    char icode;
} Residue;


/* ============================================================
   Helper function: trim whitespace
   ============================================================ */

void trim_string(char *str)
{
    char *start = str;
    char *end;

    while (*start == ' ' || *start == '\t')
        start++;

    if (start != str)
        memmove(str, start, strlen(start) + 1);

    end = str + strlen(str) - 1;

    while (end >= str &&
           (*end == ' ' || *end == '\t' ||
            *end == '\n' || *end == '\r'))
    {
        *end = '\0';
        end--;
    }
}


/* ============================================================
   Helper function: extract PDB element
   ============================================================ */

void get_element(const char *line, char *element)
{
    element[0] = '\0';
    element[1] = '\0';
    element[2] = '\0';

    /*
       PDB element berada pada kolom 77-78,
       yang pada array C adalah index 76-77.
    */

    if (strlen(line) >= 78)
    {
        element[0] = line[76];
        element[1] = line[77];
        element[2] = '\0';

        trim_string(element);
    }

    /*
       Jika element kosong, gunakan karakter pertama
       dari atom name sebagai fallback.
    */

    if (element[0] == '\0')
    {
        char atom_name[5];

        strncpy(atom_name, line + 12, 4);
        atom_name[4] = '\0';

        trim_string(atom_name);

        if (atom_name[0] != '\0')
        {
            element[0] = atom_name[0];
            element[1] = '\0';
        }
    }

    /*
       Ubah menjadi uppercase
    */

    for (int i = 0; element[i] != '\0'; i++)
    {
        if (element[i] >= 'a' && element[i] <= 'z')
        {
            element[i] =
                element[i] - 'a' + 'A';
        }
    }
}


/* ============================================================
   Check apakah atom merupakan hydrogen
   ============================================================ */

int is_hydrogen(const char *element)
{
    return element[0] == 'H';
}


/* ============================================================
   Cari residue
   ============================================================ */

int find_residue(
    Residue *residues,
    int residue_count,
    const char *res_name,
    char chain,
    int res_seq,
    char icode
)
{
    for (int i = 0; i < residue_count; i++)
    {
        if (residues[i].chain == chain &&
            residues[i].res_seq == res_seq &&
            residues[i].icode == icode)
        {
            return i;
        }
    }

    return -1;
}


/* ============================================================
   Main
   ============================================================ */

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Usage: %s <pdb_file> [output_file]\n", argv[0]);
        return 1;
    }

    const char *pdb_file = argv[1];

    const char *output_file =
        (argc >= 3)
        ? argv[2]
        : "contact_map_cpu.bin";


    /* --------------------------------------------------------
       Alokasi memory
       -------------------------------------------------------- */

    Atom *atoms =
        (Atom *)malloc(
            MAX_ATOMS * sizeof(Atom)
        );

    Residue *residues =
        (Residue *)malloc(
            MAX_RESIDUES * sizeof(Residue)
        );

    if (atoms == NULL || residues == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }


    /* --------------------------------------------------------
       Buka PDB
       -------------------------------------------------------- */

    FILE *file = fopen(pdb_file, "r");

    if (file == NULL)
    {
        printf("Tidak dapat membuka file: %s\n", pdb_file);

        free(atoms);
        free(residues);

        return 1;
    }


    /* --------------------------------------------------------
       Parsing
       -------------------------------------------------------- */

    char line[256];

    int atom_count = 0;
    int residue_count = 0;

    while (fgets(line, sizeof(line), file))
    {
        if (strncmp(line, "ATOM  ", 6) != 0)
            continue;


        /* Alternate location */
        char alt_loc = line[16];

        if (alt_loc != ' ' && alt_loc != 'A')
            continue;


        /* Element */
        char element[3];

        get_element(line, element);


        /* Ignore hydrogen */

        if (is_hydrogen(element))
            continue;


        if (atom_count >= MAX_ATOMS)
        {
            printf("MAX_ATOMS exceeded.\n");
            fclose(file);

            free(atoms);
            free(residues);

            return 1;
        }


        /* ----------------------------------------------------
           Parse fields
           ---------------------------------------------------- */

        Atom *atom = &atoms[atom_count];

        memset(atom, 0, sizeof(Atom));

        strncpy(atom->atom_name, line + 12, 4);
        atom->atom_name[4] = '\0';
        trim_string(atom->atom_name);

        strncpy(atom->res_name, line + 17, 3);
        atom->res_name[3] = '\0';
        trim_string(atom->res_name);

        atom->chain = line[21];

        char temp[16];

        memset(temp, 0, sizeof(temp));

        strncpy(temp, line + 22, 4);
        atom->res_seq = atoi(temp);

        atom->icode = line[26];

        memset(temp, 0, sizeof(temp));
        strncpy(temp, line + 30, 8);
        atom->x = atof(temp);

        memset(temp, 0, sizeof(temp));
        strncpy(temp, line + 38, 8);
        atom->y = atof(temp);

        memset(temp, 0, sizeof(temp));
        strncpy(temp, line + 46, 8);
        atom->z = atof(temp);

        strcpy(atom->element, element);


        /* ----------------------------------------------------
           Cari / buat residue
           ---------------------------------------------------- */

        int residue_index =
            find_residue(
                residues,
                residue_count,
                atom->res_name,
                atom->chain,
                atom->res_seq,
                atom->icode
            );


        if (residue_index == -1)
        {
            if (residue_count >= MAX_RESIDUES)
            {
                printf("MAX_RESIDUES exceeded.\n");
                fclose(file);

                free(atoms);
                free(residues);

                return 1;
            }

            residue_index = residue_count;

            strcpy(
                residues[residue_count].res_name,
                atom->res_name
            );

            residues[residue_count].chain =
                atom->chain;

            residues[residue_count].res_seq =
                atom->res_seq;

            residues[residue_count].icode =
                atom->icode;

            residue_count++;
        }


        atom->residue_index = residue_index;

        atom_count++;
    }

    fclose(file);


    /* --------------------------------------------------------
       Informasi dataset
       -------------------------------------------------------- */

    printf("\n");
    printf("============================================\n");
    printf("CPU CONTACT MAP\n");
    printf("============================================\n");

    printf("PDB file       : %s\n", pdb_file);
    printf("Heavy atoms    : %d\n", atom_count);
    printf("Residues       : %d\n", residue_count);


    long long total_pairs =
        ((long long)atom_count *
         (atom_count - 1)) / 2;

    printf(
        "Atom pairs     : %lld\n",
        total_pairs
    );

    printf(
        "Cutoff         : %.2f Angstrom\n",
        CUTOFF
    );


    /* --------------------------------------------------------
       Allocate contact map
       -------------------------------------------------------- */

    size_t map_size =
        (size_t)residue_count *
        residue_count;

    unsigned char *contact_map =
        (unsigned char *)calloc(
            map_size,
            sizeof(unsigned char)
        );


    if (contact_map == NULL)
    {
        printf("Contact map allocation failed.\n");

        free(atoms);
        free(residues);

        return 1;
    }


    /* --------------------------------------------------------
       Pairwise distance calculation
       -------------------------------------------------------- */

    double cutoff_squared =
        CUTOFF * CUTOFF;

    long long contact_atom_pairs = 0;


    printf("\nStarting pairwise calculation...\n");

    clock_t start = clock();


    for (int i = 0; i < atom_count; i++)
    {
        const Atom *a = &atoms[i];

        for (int j = i + 1;
             j < atom_count;
             j++)
        {
            const Atom *b = &atoms[j];


            /* Atom dari residue yang sama tidak dihitung */

            if (a->residue_index ==
                b->residue_index)
            {
                continue;
            }


            double dx =
                (double)a->x - b->x;

            double dy =
                (double)a->y - b->y;

            double dz =
                (double)a->z - b->z;


            double distance_squared =
                dx * dx +
                dy * dy +
                dz * dz;


            /* Contact */

            if (distance_squared <=
                cutoff_squared)
            {
                int r1 =
                    a->residue_index;

                int r2 =
                    b->residue_index;


                contact_map[
                    (size_t)r1 *
                    residue_count +
                    r2
                ] = 1;

                contact_map[
                    (size_t)r2 *
                    residue_count +
                    r1
                ] = 1;


                contact_atom_pairs++;
            }
        }
    }


    clock_t end = clock();


    double elapsed =
        (double)(end - start)
        / CLOCKS_PER_SEC;


    /* --------------------------------------------------------
       Hitung residue contact
       -------------------------------------------------------- */

    long long contact_residue_pairs = 0;


    for (int i = 0;
         i < residue_count;
         i++)
    {
        for (int j = i + 1;
             j < residue_count;
             j++)
        {
            if (contact_map[
                    (size_t)i *
                    residue_count +
                    j
                ])
            {
                contact_residue_pairs++;
            }
        }
    }


    /* --------------------------------------------------------
       Simpan contact map
       -------------------------------------------------------- */

    FILE *output =
        fopen(output_file, "wb");


    if (output == NULL)
    {
        printf(
            "Tidak dapat membuat output: %s\n",
            output_file
        );

        free(contact_map);
        free(atoms);
        free(residues);

        return 1;
    }


    fwrite(
        contact_map,
        sizeof(unsigned char),
        map_size,
        output
    );


    fclose(output);


    /* --------------------------------------------------------
       Hasil
       -------------------------------------------------------- */

    printf("\n");
    printf("============================================\n");
    printf("RESULT\n");
    printf("============================================\n");

    printf(
        "Atom pairs checked       : %lld\n",
        total_pairs
    );

    printf(
        "Contact atom pairs       : %lld\n",
        contact_atom_pairs
    );

    printf(
        "Contact residue pairs    : %lld\n",
        contact_residue_pairs
    );

    printf(
        "CPU execution time       : %.6f seconds\n",
        elapsed
    );

    printf(
        "Contact map output       : %s\n",
        output_file
    );


    /*--------------------------------------------------------
       Cleanup
       -------------------------------------------------------- */

    free(contact_map);
    free(atoms);
    free(residues);

    return 0;
}