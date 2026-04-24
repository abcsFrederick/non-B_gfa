PGMNAME = gfa
#CC = icc
CC = gcc -lm
OBJ = $(PGMNAME).o \
cdna.o     findIR.o    nulls.o           process_repeats.o \
findAPR.o  findMR.o    print_gff_file.o  rcdna.o \
findDR.o   findSTR.o   is_subset.o  print_tsv_file.o  read_fasta.o \
findGQ.o   findZDNA.o  print_usage.c     read_mult_fasta.o
CFLAGS = -O2 
#CFLAGS = -O2 -ansi
LFLAGS = 
$(PGMNAME): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LFLAGS)

clean:
	rm *.o
