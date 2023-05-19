n = [8, 16, 32, 64];

figure(1);
onenorm = [0.0145325,0.0148297,0.0149156,0.0149388];
twonorm = [0.00202271, 0.00107451, 0.000554058, 0.000281361];
infnorm = [0.00043333, 0.00012174, 3.2425e-05, 8.36143e-06];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Dirichlet Boundary Condition");
legend("1-norm","2-norm","inf-norm");

figure(2);
onenorm = [0.246477, 0.213941, 0.197336, 0.188962];
twonorm = [0.0299392, 0.0144453, 0.00705501, 0.0034803];
infnorm = [0.00633653, 0.00178607, 0.000475579, 0.000122805];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Pure Neumann Boundary Condition");
legend("1-norm","2-norm","inf-norm");

figure(3);
onenorm = [0.143462, 0.0889411, 0.056838, 0.0401705];
twonorm = [0.0326952, 0.0108976, 0.00339045, 0.00106188];
infnorm = [0.0228833, 0.00664539, 0.00179044, 0.00046465];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left Edge (M1)");
legend("1-norm","2-norm","inf-norm");

figure(4);
onenorm = [0.215035, 0.147904, 0.109055, 0.0890899];
twonorm = [0.0486287, 0.0174353, 0.00624647, 0.00243048];
infnorm = [0.031503, 0.00889554, 0.00234294, 0.000598197];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Right Edge (M2)");
legend("1-norm","2-norm","inf-norm");

figure(5);
onenorm = [0.072936, 0.0415, 0.0247548, 0.0172633];
twonorm = [0.0176672, 0.00560707, 0.00166916, 0.00050921];
infnorm = [0.0108471, 0.00309184, 0.000825873, 0.00021375];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Bottom Edge (M3)");
legend("1-norm","2-norm","inf-norm");

figure(6);
onenorm = [0.202032, 0.130024, 0.0985172, 0.0874886];
twonorm = [0.0491148, 0.0163039, 0.00546685, 0.00208828];
infnorm = [0.0289028, 0.00832314, 0.00223431, 0.000579677];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Top Edge (M4)");
legend("1-norm","2-norm","inf-norm");

figure(7);
onenorm = [0.254014, 0.154123, 0.106667, 0.0879019];
twonorm = [0.051953, 0.0171458, 0.00564736, 0.00209888];
infnorm = [0.0288827, 0.00832267, 0.00223442, 0.000579703];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Bottom and Top Edge (M5)");
legend("1-norm","2-norm","inf-norm");

figure(8);
onenorm = [0.29149, 0.19015, 0.130075, 0.0991921];
twonorm = [0.0563724, 0.0197859, 0.00681676, 0.00253692];
infnorm = [0.0313607, 0.00888337, 0.00234201, 0.000598128];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left and Right Edge (M6)");
legend("1-norm","2-norm","inf-norm");

figure(9);
onenorm = [0.29212, 0.152903, 0.104804, 0.0904629];
twonorm = [0.0482946, 0.0147852, 0.00503098, 0.00208664];
infnorm = [0.0271911, 0.00727603, 0.00187866, 0.000476681];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Right and Top Edge (M7)");
legend("1-norm","2-norm","inf-norm");

figure(10);
onenorm = [0.331551, 0.173204, 0.0897964, 0.0477838];
twonorm = [0.0462259, 0.0139371, 0.00400939, 0.00114999];
infnorm = [0.0235183, 0.00669039, 0.00179303, 0.000464747];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left and Bottom Edge (M8)");
legend("1-norm","2-norm","inf-norm");

figure(11);
onenorm = [0.162944, 0.113716, 0.145256, 0.170709];
twonorm = [0.0386866, 0.0129444, 0.00627946, 0.00344622];
infnorm = [0.0296346, 0.00834688, 0.00223158, 0.000578921];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left and Top Edge (M9)");
legend("1-norm","2-norm","inf-norm");

figure(12);
onenorm = [0.644572, 0.444562, 0.338251, 0.283986];
twonorm = [0.0863696, 0.0324974, 0.0129666, 0.00560963];
infnorm = [0.03286, 0.00904178, 0.00235896, 0.000600025];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Right and Bottom Edge (M10)");
legend("1-norm","2-norm","inf-norm");

figure(13);
onenorm = [0.398261, 0.344147, 0.323138, 0.315933];
twonorm = [0.0541188, 0.0244664, 0.0119859, 0.00603356];
infnorm = [0.0263621, 0.00716928, 0.00186524, 0.000475006];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Right, Bottom, and Top Edge (M11)");
legend("1-norm","2-norm","inf-norm");

figure(14);
onenorm = [0.222654, 0.23478, 0.254294, 0.268824];
twonorm = [0.0402395, 0.0173261, 0.00894785, 0.00475327];
infnorm = [0.0292519, 0.00830343, 0.00222661, 0.000578336];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left, Bottom, and Top Edge (M12)");
legend("1-norm","2-norm","inf-norm");

figure(15);
onenorm = [0.18399, 0.149266, 0.137503, 0.135201];
twonorm = [0.0256505, 0.0113735, 0.00542243, 0.00268511];
infnorm = [0.0100771, 0.0029309, 0.00079886, 0.000209653];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left, Right, and Top Edge (M13)");
legend("1-norm","2-norm","inf-norm");

figure(16);
onenorm = [0.657373, 0.514191, 0.439643, 0.403859];
twonorm = [0.0857476, 0.035418, 0.0154803, 0.00719951];
infnorm = [0.0326801, 0.00905016, 0.00236278, 0.00060071];
hold on
plot(n,onenorm, "r--");
plot(n,twonorm, "g:");
plot(n,infnorm, "b.-");
hold off
title("Errors for Neumann BC on Left, Right, and Bottom Edge (M14)");
legend("1-norm","2-norm","inf-norm");