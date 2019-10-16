Please prepare the data file contains the following information in this order:
• Number of constraints
• Number of integer variables
• First objective function coefficients
• Second objective function coefficients
• Third objective function coefficients
• Constraints coefficient matrix for integer variables
• Right hand side value of each constraint
• A binary value indicating the type of constraint. A zero means an = constraint and a one means a ≤ constraint.

To run this program, please set argv[1] as the input data file, argv[2] as the reorte file, argv[3] as the nondominated points set file, and argv[4] as the name of instance. 
For example, .\quadrant_shrinking_method "Input_data.txt" "Report.txt" "ND_set.txt" "AP_instance".

The report file will present the instance name, run time, and number of nondominated points.
The nondominated points set file will present all of the nondominated points.
