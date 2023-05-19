%{
dir32res = [8.63615; 0.273882; 0.0101765; 0.000413073; 1.77323e-05; 7.92344e-07; 3.65258e-08; ...
    1.72759e-09; 8.41283e-11; 4.54747e-12; 9.09495e-13; 9.09495e-13; 9.09495e-13; ...
    9.09495e-13; 1.36424e-12; 9.09495e-13; 9.09495e-13; 6.82121e-13; 9.09495e-13; ...
    9.09495e-13; 6.82121e-13; ];

dir32rate = [0.968287; 0.962843; 0.959409; 0.957072; 0.955316; ...
    0.953902; 0.952702; 0.951303; 0.945946; 0.8; 0; 0; 0; ...
    -0.5; 0.333333; 0; 0.25; -0.333333; 0; 0.25; ];

dir64res = [34.5352; 1.09543; 0.0407165; 0.00165379; 7.10696e-05; ...
    3.1807e-06; 1.46945e-07; 6.96855e-09; 3.40151e-10; ...
    2.00089e-11; 3.63798e-12; 5.45697e-12; 2.72848e-12; ...
    2.72848e-12; 1.81899e-12; 1.81899e-12; 1.81899e-12; ...
    1.81899e-12; 1.81899e-12; 1.81899e-12; 1.81899e-12;];

dir64rate = [0.968281; 0.962831; 0.959383; 0.957026; 0.955245; ...
    0.953801; 0.952577; 0.951188; 0.941176; 0.818182; -0.5; 0.5; 0; ...
    0.333333; 0; 0; 0; 0; 0; 0; ];

dir128res = [138.128; 4.38121; 0.162844; 0.00661418; 0.000284235; ...
    1.27209e-05; 5.87694e-07; 2.78669e-08; 1.34605e-09; ...
    7.27596e-11; 1.45519e-11; 1.45519e-11; 1.45519e-11; ...
    2.18279e-11; 1.45519e-11; 1.45519e-11; 2.18279e-11; ...
    1.45519e-11; 1.45519e-11; 1.45519e-11; 1.45519e-11;];

dir128rate = [0.968281; 0.962831; 0.959383; 0.957026; 0.955245; ...
    0.953801; 0.952583; 0.951697; 0.945946; 0.8; 0; 0; -0.5; ...
    0.333333; 0; -0.5; 0.333333; 0; 0; 0; ];

dir256res = [552.497; 17.5243; 0.651352; 0.0264556; 0.00113688; ...
    5.08809e-05; 2.35069e-06; 1.11497e-07; 5.44242e-09; ...
    3.49246e-10; 5.82077e-11; 8.73115e-11; 5.82077e-11; ...
    8.73115e-11; 5.82077e-11; 8.73115e-11; 5.82077e-11; ...
    8.73115e-11; 5.82077e-11; 8.73115e-11; 5.82077e-11;];

dir256rate = [0.968282; 0.962831; 0.959384; 0.957027; 0.955245; ...
    0.9538; 0.952568; 0.951188; 0.935829; 0.833333; -0.5; 0.333333; -0.5; ...
    0.333333; -0.5; 0.333333; -0.5; 0.333333; -0.5; 0.333333; ];

figure(1)
plot(1:21, dir32res);
hold on
plot(1:21, dir64res);
plot(1:21, dir128res);
plot(1:21, dir256res);
hold off
title("Residual for each V-cycle of Dirichlet Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");

figure(2)
plot(2:21, dir32rate);
hold on
plot(2:21, dir64rate);
plot(2:21, dir128rate);
plot(2:21, dir256rate);
hold off
title("Reduction Rate for Each V-cycle of Dirichlet Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");
%}

%{
neu32res = [0.215111; 0.120573; 0.0264598; 0.0121078; 0.00330463; ...
    0.00144024; 0.000433704; 0.000180426; 5.74253e-05; ...
    2.3009e-05; 7.58376e-06; 2.9562e-06; 9.97479e-07; ...
    3.8131e-07; 1.30852e-07; 4.93058e-08; 1.71327e-08; ...
    6.38651e-09; 2.24001e-09; 8.28269e-10; 2.92553e-10;];

neu32rate = [0.439484; 0.78055; 0.542409; 0.727065; 0.564175; ...
    0.698867; 0.583989; 0.681723; 0.599324; 0.6704; 0.610193; 0.662581; ...
    0.617726;0.656836; 0.623194; 0.652521; 0.627234; 0.649259; 0.630239; ...
    0.64679; ];

neu64res = [0.836914; 0.558029; 0.293898; 0.176624; 0.105768; ...
    0.0665875; 0.0423169; 0.0274371; 0.0179042; 0.0117749; 0.00776875; ...
    0.00514103; 0.00340697; 0.00226043; 0.00150063; 0.000996681; 0.000662134; ...
    0.00043996; 0.000292364; 0.000194297; 0.000129129;];

neu64rate = [0.333231; 0.473328; 0.399028; 0.401171; 0.370437; ...
    0.364492; 0.351628; 0.347446; 0.342336; 0.340229; 0.338242; 0.337299; ...
    0.336528; 0.336127; 0.335827; 0.335661; 0.335543; ...
    0.335475; 0.335429; 0.335401; ];

neu128res = [1.96927; 1.83796; 1.6361; 1.51632; 1.43755; 1.39676; 1.38272; ...
    1.39025; 1.41444; 1.45208; 1.50069; 1.55854; 1.62437; 1.6973; 1.77676; ...
    1.86236; 1.95391; 2.0513; 2.15456; 2.26377; 2.37907;];

neu128rate = [0.0666839; 0.109825; 0.0732127; 0.0519492; 0.0283756; ...
    0.0100483; -0.0054443; -0.0174003; -0.026611; -0.0334786; -0.0385474; ...
    -0.042236; -0.0449009; -0.0468135; -0.0481804; -0.0491542; -0.0498464; ...
    -0.0503376; -0.0506859; -0.0509326; ];

neu256res = [4.36269; 5.34227; 7.00904; 9.40185; 12.8479; 17.8723; 25.2599; ...
    36.2011; 52.4968; 76.8758; 113.474; 168.563; 251.656; 377.185; 567.052; ...
    854.493; 1289.96; 1950.01; 2950.9; 4469.08; 6772.41;];

neu256rate = [-0.224535; -0.311998; -0.341389; -0.366532; -0.391063; ...
    -0.413356; -0.433145; -0.450145; -0.46439; ...
    -0.476069; -0.485477; -0.492948; -0.498813; ...
    -0.503378; -0.506904; -0.509615; -0.51169; ...
    -0.513273; -0.514478; -0.515393; ];

figure(1)
plot(1:21, neu32res);
hold on
plot(1:21, neu64res);
plot(1:21, neu128res);
plot(1:21, neu256res);
hold off
title("Residual for Each V-cycle of Pure Neumann Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");

figure(2)
plot(2:21, neu32rate);
hold on
plot(2:21, neu64rate);
plot(2:21, neu128rate);
plot(2:21, neu256rate);
hold off
title("Reduction Rate for Each V-cycle of Pure Neumann Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");
%}

%{
left32res = [0.185098; 0.0959573; 0.0155313; 0.00778045; 0.00166602; ...
    0.000839214; 0.000213487; 0.000100915; 2.8499e-05; ...
    1.26133e-05; 3.81285e-06; 1.60107e-06; 5.07049e-07; ...
    2.04805e-07; 6.70377e-08; 2.63219e-08; 8.82782e-09; ...
    3.39395e-09; 1.15911e-09; 4.3865e-10; 1.51864e-10;];

left32rate = [0.481587; 0.838143; 0.499048; 0.785871; 0.496277; ...
    0.745611; 0.527301; 0.717594; 0.557412; ...
    0.697712; 0.580085; 0.683306; 0.596084; ...
    0.672676; 0.607356; 0.664621; 0.61554; 0.658478; 0.621562; 0.653793; ];

left64res = [0.835107; 0.531492; 0.265669; 0.153434; 0.0884876; ...
    0.0543205; 0.0338342; 0.021661; 0.014007; ...
    0.00916032; 0.0060207; 0.00397489; 0.00263008; ...
    0.00174331; 0.00115662; 0.000767896; 0.000510017; ...
    0.000338831; 0.000225139; 0.000149612; 9.94279e-05; ];

left64rate = [0.363564; 0.500146; 0.422459; 0.423287; 0.386123; ...
    0.377138; 0.359788; 0.353355; 0.346019; ...
    0.342741; 0.339796; 0.338326; 0.337165; ...
    0.33654; 0.336084; 0.335826; 0.335647; 0.335541; 0.33547; 0.335428; ];

left128res = [1.96613; 1.80713; 1.58055; 1.44082; 1.34545; ...
    1.29046; 1.26392; 1.26011; 1.27371; 1.30121; 1.33991; 1.38787; 1.44372; ...
    1.50646; 1.57543; 1.65018; 1.73042; 1.81603; 1.90697; 2.00327; 2.10503;];

left128rate = [0.0808656; 0.125381; 0.0884081; 0.0661947; 0.0408683; ...
    0.0205641; 0.00301893; -0.0107948; -0.0215921; ...
    -0.0297369; -0.0357993; -0.0402389; -0.043461; ...
    -0.0457811; -0.0474432; -0.0486293; -0.0494734; ...
    -0.050073; -0.0504984; -0.0507998; ];

left256res = [4.36853; 5.32168; 6.93611; 9.2429; 12.5513; ...
    17.3585; 24.4071; 34.8229; 50.3088; 73.4441; 108.138; 160.317; 238.971; ...
    357.737; 537.307; 809.083; 1220.73; 1844.58; 2790.45; 4225.03; 6401.39;];

left256rate = [-0.218185; -0.303369; -0.332577; -0.35794; -0.383004; ...
    -0.406061; -0.426754; -0.444702; -0.459866; -0.472385; -0.482525; -0.490614; -0.496988; ...
    -0.501961; -0.505812; -0.508777; -0.511049; -0.512785; -0.514106; -0.515111;];

figure(1)
plot(1:21, left32res);
hold on
plot(1:21, left64res);
plot(1:21, left128res);
plot(1:21, left256res);
hold off
title("Residual for Each V-cycle of Left Neumann Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");

figure(2)
plot(2:21, left32rate);
hold on
plot(2:21, left64rate);
plot(2:21, left128rate);
plot(2:21, left256rate);
hold off
title("Reduction Rate for Each V-cycle of Left Neumann Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");
%}

right32res = [0.214438; 0.120358; 0.0263927; 0.0120853; 0.00329679; ...
    0.00143749; 0.000432714; 0.000180073; 5.72983e-05; ...
    2.29631e-05; 7.56721e-06; 2.95024e-06; 9.95323e-07; ...
    3.80534e-07; 1.30571e-07; 4.92047e-08; 1.70962e-08; ...
    6.37334e-09; 2.23526e-09; 8.26549e-10; 2.91934e-10;];

right32rate = [0.438727; 0.780715; 0.542095; 0.727208; 0.563973; ...
    0.698979; 0.583854; 0.681804; 0.599235; 0.670463; 0.610128; 0.66263; 0.617678; ...
    0.656873; 0.623159; 0.652549; 0.627207; 0.64928; 0.630222; 0.646805; ];

right64res = [0.836497; 0.557779; 0.293757; 0.17654; 0.105716; ...
    0.0665547; 0.0422958; 0.0274234; 0.0178952; ...
    0.011769; 0.00776484; 0.00513844; 0.00340525; ...
    0.00225929; 0.00149988; 0.000996179; 0.0006618; ...
    0.000439738; 0.000292217; 0.000194199; 0.000129064;];

right64rate = [0.333197; 0.473346; 0.399026; 0.401179; 0.370438; ...
    0.364496; 0.351629; 0.347448; 0.342336; 0.34023; 0.338242; 0.337299; 0.336528; ...
    0.336127; 0.335827; 0.335661; 0.335543; 0.335475; 0.335429; 0.335401; ];

right128res = [1.96906; 1.83776; 1.63592; 1.51615; 1.43738; ...
    1.39659; 1.38255; 1.39008; 1.41427; 1.4519; 1.50051; 1.55835; 1.62416; ...
    1.69709; 1.77654; 1.86213; 1.95366; 2.05104; 2.15429; 2.26348; 2.37877;];

right128rate = [0.066684; 0.109828; 0.0732152; 0.0519519; 0.0283777; ...
    0.0100501; -0.00544287; -0.0173992; -0.0266102; ...
    -0.0334779; -0.038547; -0.0422357; -0.0449007; ...
    -0.0468133; -0.0481803; -0.0491541; -0.0498463; -0.0503376; -0.0506859; -0.0509325; ];

right256res = [4.36258; 5.34213; 7.00886; 9.4016; 12.8476; ...
    17.8718; 25.2592; 36.2001; 52.4953; 76.8736; 113.471; 168.558; 251.648; ...
    377.174; 567.035; 854.467; 1289.92; 1949.95; 2950.81; 4468.94; 6772.2;];

right256rate = [-0.224535; -0.311997; -0.341388; -0.366531; -0.391062; ...
    -0.413355; -0.433145; -0.450144; -0.464389; -0.476068; -0.485477; ...
    -0.492948; -0.498813; -0.503377; -0.506904; -0.509615; -0.51169; ...
    -0.513273; -0.514478; -0.515393; ];

figure(1)
plot(1:21, left32res);
hold on
plot(1:21, left64res);
plot(1:21, left128res);
plot(1:21, left256res);
hold off
title("Residual for Each V-cycle of Right Neumann Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");

figure(2)
plot(2:21, left32rate);
hold on
plot(2:21, left64rate);
plot(2:21, left128rate);
plot(2:21, left256rate);
hold off
title("Reduction Rate for Each V-cycle of Right Neumann Boundary Condition");
legend("n = 32", "n = 64", "n = 128", "n = 256");