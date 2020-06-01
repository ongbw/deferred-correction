function [A,b,c,s,type,order] = init_integrator(integrator)

% type:
% 1 - explicit
% 2 - diagonally implicit
% 3 - fully implicit

switch integrator
  case{1}
    % Forward Euler
    type = 0; % explicit
    s = 1; % number of stages
    A = [0];
    c = 0;
    b = 1;
    order = 1;
    
  case{2}
    % Midpoint
    type = 0; % explicit
    s = 2;
    A = [0 0 ; 
               1/2 0];
    c = [0; 0.5];
    b = [0 1];
    order = 2;
    
  case{3}
    % Trapezoid 
    type = 1; % implicit
    s = 2;
    A = [0 0;
               1/2 1/2];
    c = [0; 1];
    b = [1/2 1/2];
    order = 2;
    
  case{4}
    % RK3
    type = 0; % explicit
    s = 3;
    A = [0 0 0;
               0.5 0 0;
               -1 2 0];
    c = [0;0.5; 1];
    b = [1/6 2/3 1/6];
    order = 3;
    
  case{5}
    % RK4
    type = 0; % explicit
    s = 4;
    A = [0 0 0 0;
               0.5 0 0 0;
               0 0.5 0 0;
               0 0 1 0];
    c = [0; 0.5; 0.5; 1];
    b = [1/6 1/3 1/3 1/6];
    order = 4;
    
  case{6}
    % Backward Euler
    type = 1; % implicit
    s = 1;
    A = 1;
    c = 1;
    b = 1;
    order = 1;
    
  case{7}
    % SDIRK 2
    type = 1; % implicit
    gamma = (2-sqrt(2))/2;
    s = 2;
    A = [gamma 0;
               1-gamma gamma];
    c = [gamma; 1];
    b = [(1-gamma) gamma];
    order = 2;
    
  case{8}
    % built in ode45 integrator
    order = 4;
    type  = 2; % built in
    A = [];
    s = [];
    c = [];
    b = [];
    
  case{9}
    % RK5
    order = 5;
    type = 0;
    s = 7;
   
    c = [0;1/5;3/10;4/5;8/9;1;1];
    A = zeros(s,s);
    A(2,1) = 1/5;
    A(3,1) = 3/40;
    A(3,2) = 9/40;
    A(4,1) = 44/45;
    A(4,2) = -56/15;
    A(4,3) = 32/9;
    A(5,1) = 19372/6561;
    A(5,2) = -25360/2187;
    A(5,3) = 64448/6561;
    A(5,4) = -212/729;
    A(6,1) = 9017/3168;
    A(6,2) = -355/33;
    A(6,3) = 46732/5247;
    A(6,4) = 49/176;
    A(6,5) = -5103/18656;
    A(7,1) = 35/384;
    A(7,2) = 0;
    A(7,3) = 500/1113;
    A(7,4) = 125/192;
    A(7,5) = -2187/6784;
    A(7,6) = 11/84;
    
    b = A(7,:);
    
  case{10}
    % RK6
    order = 6;
    type = 0;
    s = 10;
    A = zeros(s,s);
    c = zeros(s,1);
    c(1) =  0;
    c(2) = .5e-2;
    c(3) = .1088888888888888888888888888888888888889;
    c(4) = .1633333333333333333333333333333333333333;
    c(5) = .4555;
    c(6) = .6095094489978381317087004421486024949638;
    c(7) = .884;
    c(8) = .925;
    c(9) = 1.;
    c(10) = 1.;
    
    A(2,1) =  0.5e-2;

    A(3,1) = -1.076790123456790123456790123456790123457;
    A(3,2) =  1.185679012345679012345679012345679012346;
    
    A(4,1) =  .4083333333333333333333333333333333333333e-1;
    A(4,2) =  0.;
    A(4,3) =  .1225;
    
    A(5,1) =  .6389139236255726780508121615993336109954;
    A(5,2) =  0.;
    A(5,3) = -2.455672638223656809662640566430653894211;
    A(5,4) =  2.272258714598084131611828404831320283215;
    
    A(6,1) = -2.661577375018757131119259297861818119279;
    A(6,2) =  0.;
    A(6,3) =  10.80451388645613769565396655365532838482;
    A(6,4) = -8.353914657396199411968048547819291691541;
    A(6,5) =  .8204875949566569791420417341743839209619;

    A(7,1) =  6.067741434696770992718360183877276714679;
    A(7,2) =  0.;
    A(7,3) = -24.71127363591108579734203485290746001803;
    A(7,4) =  20.42751793078889394045773111748346612697;
    A(7,5) = -1.906157978816647150624096784352757010879;
    A(7,6) =  1.006172249242068014790040335899474187268;

    A(8,1) =  12.05467007625320299509109452892778311648;
    A(8,2) =  0.;
    A(8,3) = -49.75478495046898932807257615331444758322;
    A(8,4) =  41.14288863860467663259698416710157354209;
    A(8,5) = -4.461760149974004185641911603484815375051;
    A(8,6) =  2.042334822239174959821717077708608543738;
    A(8,7) = -0.9834843665406107379530801693870224403537e-1;

    A(9,1) =  10.13814652288180787641845141981689030769;
    A(9,2) =  0.;
    A(9,3) = -42.64113603171750214622846006736635730625;
    A(9,4) =  35.76384003992257007135021178023160054034;
    A(9,5) = -4.348022840392907653340370296908245943710;
    A(9,6) =  2.009862268377035895441943593011827554771;;
    A(9,7) =  .3487490460338272405953822853053145879140;
    A(9,8) = -.2714390051048312842371587140910297407572;

    A(10,1) = -45.03007203429867712435322405073769635151;
    A(10,2) =  0.;
    A(10,3) =  187.3272437654588840752418206154201997384;
    A(10,4) = -154.0288236935018690596728621034510402582;
    A(10,5) =  18.56465306347536233859492332958439136765;
    A(10,6) = -7.141809679295078854925420496823551192821;
    A(10,7) =  1.308808578161378625114762706007696696508;
    A(10,8) =  0.;
    A(10,9) =  0.;

    % sixth order weights
    b(1) =  .4460860660634117628731817597479197781432e-1;
    b(2) =  0.;
    b(3) =  0.;
    b(4) =  .2671640378571372680509102260943837899738;
    b(5) =  .2201018300177293019979715776650753096323;
    b(6) =  .2188431703143156830983120833512893824578;
    b(7) =  .2289871705411202883378173889763552365362;
    b(8) =  0.;
    b(9) =  0.;
    b(10) =  .2029518466335628222767054793810430358554e-1;

  case{11}
    % RK 7 
    order = 7;
    type = 0;
    s = 10;
    A = zeros(s,s);
    c = zeros(s,1);
    c(1) =  0;
    c(2) = .5e-2;
    c(3) = .1088888888888888888888888888888888888889;
    c(4) = .1633333333333333333333333333333333333333;
    c(5) = .4555;
    c(6) = .6095094489978381317087004421486024949638;
    c(7) = .884;
    c(8) = .925;
    c(9) = 1.;
    c(10) = 1.;
    
    A(2,1) =  0.5e-2;

    A(3,1) = -1.076790123456790123456790123456790123457;
    A(3,2) =  1.185679012345679012345679012345679012346;
    
    A(4,1) =  .4083333333333333333333333333333333333333e-1;
    A(4,2) =  0.;
    A(4,3) =  .1225;
    
    A(5,1) =  .6389139236255726780508121615993336109954;
    A(5,2) =  0.;
    A(5,3) = -2.455672638223656809662640566430653894211;
    A(5,4) =  2.272258714598084131611828404831320283215;
    
    A(6,1) = -2.661577375018757131119259297861818119279;
    A(6,2) =  0.;
    A(6,3) =  10.80451388645613769565396655365532838482;
    A(6,4) = -8.353914657396199411968048547819291691541;
    A(6,5) =  .8204875949566569791420417341743839209619;

    A(7,1) =  6.067741434696770992718360183877276714679;
    A(7,2) =  0.;
    A(7,3) = -24.71127363591108579734203485290746001803;
    A(7,4) =  20.42751793078889394045773111748346612697;
    A(7,5) = -1.906157978816647150624096784352757010879;
    A(7,6) =  1.006172249242068014790040335899474187268;

    A(8,1) =  12.05467007625320299509109452892778311648;
    A(8,2) =  0.;
    A(8,3) = -49.75478495046898932807257615331444758322;
    A(8,4) =  41.14288863860467663259698416710157354209;
    A(8,5) = -4.461760149974004185641911603484815375051;
    A(8,6) =  2.042334822239174959821717077708608543738;
    A(8,7) = -0.9834843665406107379530801693870224403537e-1;

    A(9,1) =  10.13814652288180787641845141981689030769;
    A(9,2) =  0.;
    A(9,3) = -42.64113603171750214622846006736635730625;
    A(9,4) =  35.76384003992257007135021178023160054034;
    A(9,5) = -4.348022840392907653340370296908245943710;
    A(9,6) =  2.009862268377035895441943593011827554771;;
    A(9,7) =  .3487490460338272405953822853053145879140;
    A(9,8) = -.2714390051048312842371587140910297407572;

    A(10,1) = -45.03007203429867712435322405073769635151;
    A(10,2) =  0.;
    A(10,3) =  187.3272437654588840752418206154201997384;
    A(10,4) = -154.0288236935018690596728621034510402582;
    A(10,5) =  18.56465306347536233859492332958439136765;
    A(10,6) = -7.141809679295078854925420496823551192821;
    A(10,7) =  1.308808578161378625114762706007696696508;
    A(10,8) =  0.;
    A(10,9) =  0.;

    % seventh order weights
    b(1) =  .4715561848627222170431765108838175679569e-1;
    b(2) =  0.;
    b(3) =  0.;
    b(4) =  .2575056429843415189596436101037687580986;
    b(5) =  .2621665397741262047713863095764527711129;
    b(6) =  .1521609265673855740323133199165117535523;
    b(7) =  .4939969170032484246907175893227876844296;
    b(8) = -.2943031171403250441557244744092703429139;
    b(9) =  .8131747232495109999734599440136761892478e-1;
    b(10) =  0.;

  case{12}
    % RK 8 
    order = 8;
    type = 0;
    s = 16;
    A = zeros(s,s);
    c = zeros(s,1);
    
    c(1) = 0.;
    c(2) = .3462e-1;
    c(3) = .9702435063878044594828361677100617517633e-1;
    c(4) = .1455365259581706689224254251565092627645;
    c(5) = .561;
    c(6) = .2290079115904850126662751771814700052182;
    c(7) = .5449920884095149873337248228185299947818;
    c(8) = .645;
    c(9) = .4837500000000000000000000000000000000000;
    c(10) = .6757e-1;
    c(11) = .2500;
    c(12) = .6590650618730998549405331618649220295334;
    c(13) = .8206;
    c(14) = .9012;
    c(15) = 1.;
    c(16) = 1.;

    A(2,1) =  .3462e-1;

    A(3,1) = -.389335438857287327017042687229284478532e-1;
    A(3,2) =  .1359578945245091786499878854939346230295;

    A(4,1) =  .3638413148954266723060635628912731569111e-1;
    A(4,2) =  0.;
    A(4,3) =  .1091523944686280016918190688673819470733;

    A(5,1) =  2.025763914393969636805657604282571047511;
    A(5,2) =  0.;
    A(5,3) = -7.638023836496292020387602153091964592952;
    A(5,4) =  6.173259922102322383581944548809393545442;

    A(6,1) =  .5112275589406060872792270881648288397197e-1;
    A(6,2) =  0.;
    A(6,3) =  0.;
    A(6,4) =  .1770823794555021537929910813839068684087;
    A(6,5) =  .80277624092225014536138698108025283759e-3;

    A(7,1) =  .1316006357975216279279871693164256985334;
    A(7,2) =  0.;
    A(7,3) =  0.;
    A(7,4) = -.2957276252669636417685183174672273730699;
    A(7,5) =  .878137803564295237421124704053886667082e-1;
    A(7,6) =  .6213052975225274774321435005639430026100;

    A(8,1) =  .7166666666666666666666666666666666666667e-1;
    A(8,2) =  0.;
    A(8,3) =  0.;
    A(8,4) =  0.;
    A(8,5) =  0.;
    A(8,6) =  .3305533578915319409260346730051472207728;
    A(8,7) =  .2427799754418013924072986603281861125606;

    A(9,1) =  .7180664062500000000000000000000000000000e-1;
    A(9,2) =  0.;
    A(9,3) =  0.;
    A(9,4) =  0.;
    A(9,5) =  0.;
    A(9,6) =  .3294380283228177160744825466257672816401;
    A(9,7) =  .1165190029271822839255174533742327183599;
    A(9,8) = -.3401367187500000000000000000000000000000e-1;

    A(10,1) =  .4836757646340646986611287718844085773549e-1;
    A(10,2) =  0.;
    A(10,3) =  0.;
    A(10,4) =  0.;
    A(10,5) =  0.;
    A(10,6) =  .3928989925676163974333190042057047002852e-1;
    A(10,7) =  .1054740945890344608263649267140088017604;
    A(10,8) = -.2143865284648312665982642293830533996214e-1;
    A(10,9) = -.1041229174627194437759832813847147895623;

    A(11,1) = -.2664561487201478635337289243849737340534e-1;
    A(11,2) =  0.;
    A(11,3) =  0.;
    A(11,4) =  0.;
    A(11,5) =  0.;
    A(11,6) =  .3333333333333333333333333333333333333333e-1;
    A(11,7) = -.1631072244872467239162704487554706387141;
    A(11,8) =  .3396081684127761199487954930015522928244e-1;
    A(11,9) =  .1572319413814626097110769806810024118077;
    A(11,10) =  .2152267478031879552303534778794770376960;

    A(12,1) =  .3689009248708622334786359863227633989718e-1;
    A(12,2) =  0.;
    A(12,3) =  0.;
    A(12,4) =  0.;
    A(12,5) =  0.;
    A(12,6) = -.1465181576725542928653609891758501156785;
    A(12,7) =  .2242577768172024345345469822625833796001;
    A(12,8) =  .2294405717066072637090897902753790803034e-1;
    A(12,9) = -.35850052905728761357394424889330334334e-2;
    A(12,10) =  .8669223316444385506869203619044453906053e-1;
    A(12,11) =  .4383840651968337846196219974168630120572;

    A(13,1) = -.4866012215113340846662212357570395295088;
    A(13,2) =  0.;
    A(13,3) =  0.;
    A(13,4) =  0.;;
    A(13,5) =  0.;
    A(13,6) = -6.304602650282852990657772792012007122988;
    A(13,7) = -.281245618289472564778284183790118418111;
    A(13,8) = -2.679019236219849057687906597489223155566;
    A(13,9) =  .518815663924157511565311164615012522024;
    A(13,10) =  1.365353187603341710683633635235238678626;
    A(13,11) =  5.885091088503946585721274891680604830712;
    A(13,12) =  2.802808786272062889819965117517532194812;

    A(14,1) =  .4185367457753471441471025246471931649633;
    A(14,2) =  0.;
    A(14,3) =  0.;
    A(14,4) =  0.;
    A(14,5) =  0.;
    A(14,6) =  6.724547581906459363100870806514855026676;
    A(14,7) = -.425444280164611790606983409697113064616;
    A(14,8) =  3.343279153001265577811816947557982637749;
    A(14,9) =  .617081663117537759528421117507709784737;
    A(14,10) = -.929966123939932833937749523988800852013;
    A(14,11) = -6.099948804751010722472962837945508844846;
    A(14,12) = -3.002206187889399044804158084895173690015;
    A(14,13) =  .2553202529443445472336424602988558373637;

    A(15,1) = -.779374086122884664644623040843840506343;
    A(15,2) =  0.;
    A(15,3) =  0.;
    A(15,4) =  0.;
    A(15,5) =  0.;
    A(15,6) = -13.93734253810777678786523664804936051203;
    A(15,7) =  1.252048853379357320949735183924200895136;
    A(15,8) = -14.69150040801686878191527989293072091588;
    A(15,9) = -.494705058533141685655191992136962873577;
    A(15,10) =  2.242974909146236657906984549543692874755;
    A(15,11) =  13.36789380382864375813864978592679139881;
    A(15,12) =  14.39665048665068644512236935340272139005;
    A(15,13) = -.7975813331776800379127866056663258667437;
    A(15,14) =  .4409353709534277758753793068298041158235;

    A(16,1) =  2.058051337466886442151242368989994043993;
    A(16,2) =  0.;
    A(16,3) =  0.;
    A(16,4) =  0.;
    A(16,5) =  0.;
    A(16,6) =  22.35793772796803295519317565842520212899;
    A(16,7) =  .90949810997556332745009198137971890783;
    A(16,8) =  35.89110098240264104710550686568482456493;
    A(16,9) = -3.442515027624453437985000403608480262211;
    A(16,10) = -4.865481358036368826566013387928704014496;
    A(16,11) = -18.90980381354342625688427480879773032857;
    A(16,12) = -34.26354448030451782929251177395134170515;
    A(16,13) =  1.264756521695642578827783499806516664686;
    A(16,14) =  0.;
    A(16,15) =  0.;

    b(1) =  .1996996514886773085518508418098868756464e-1;
    b(2) =  0.;
    b(3) =  0.;
    b(4) =  0.;
    b(5) =  0.;
    b(6) =  0.;;
    b(7) =  0.;
    b(8) =  2.191499304949330054530747099310837524864;
    b(9) =  .8857071848208438030833722031786358862953e-1;
    b(10) =  .1140560234865965622484956605091432032674;
    b(11) =  .2533163805345107065564577734569651977347;
    b(12) = -2.056564386240941011158999594595981300493;
    b(13) =  .3408096799013119935160094894224543812830;
    b(14) =  0.;
    b(15) =  0.;
    b(16) =  .4834231373823958314376726739772871714902e-1;
    
 case{13}
    % RK 9 
    order = 9;
    type = 0;
    s = 16;
    A = zeros(s,s);
    c = zeros(s,1);
    
    c(1) = 0.;
    c(2) = .3462e-1;
    c(3) = .9702435063878044594828361677100617517633e-1;
    c(4) = .1455365259581706689224254251565092627645;
    c(5) = .561;
    c(6) = .2290079115904850126662751771814700052182;
    c(7) = .5449920884095149873337248228185299947818;
    c(8) = .645;
    c(9) = .4837500000000000000000000000000000000000;
    c(10) = .6757e-1;
    c(11) = .2500;
    c(12) = .6590650618730998549405331618649220295334;
    c(13) = .8206;
    c(14) = .9012;
    c(15) = 1.;
    c(16) = 1.;

    A(2,1) =  .3462e-1;

    A(3,1) = -.389335438857287327017042687229284478532e-1;
    A(3,2) =  .1359578945245091786499878854939346230295;

    A(4,1) =  .3638413148954266723060635628912731569111e-1;
    A(4,2) =  0.;
    A(4,3) =  .1091523944686280016918190688673819470733;

    A(5,1) =  2.025763914393969636805657604282571047511;
    A(5,2) =  0.;
    A(5,3) = -7.638023836496292020387602153091964592952;
    A(5,4) =  6.173259922102322383581944548809393545442;

    A(6,1) =  .5112275589406060872792270881648288397197e-1;
    A(6,2) =  0.;
    A(6,3) =  0.;
    A(6,4) =  .1770823794555021537929910813839068684087;
    A(6,5) =  .80277624092225014536138698108025283759e-3;

    A(7,1) =  .1316006357975216279279871693164256985334;
    A(7,2) =  0.;
    A(7,3) =  0.;
    A(7,4) = -.2957276252669636417685183174672273730699;
    A(7,5) =  .878137803564295237421124704053886667082e-1;
    A(7,6) =  .6213052975225274774321435005639430026100;

    A(8,1) =  .7166666666666666666666666666666666666667e-1;
    A(8,2) =  0.;
    A(8,3) =  0.;
    A(8,4) =  0.;
    A(8,5) =  0.;
    A(8,6) =  .3305533578915319409260346730051472207728;
    A(8,7) =  .2427799754418013924072986603281861125606;

    A(9,1) =  .7180664062500000000000000000000000000000e-1;
    A(9,2) =  0.;
    A(9,3) =  0.;
    A(9,4) =  0.;
    A(9,5) =  0.;
    A(9,6) =  .3294380283228177160744825466257672816401;
    A(9,7) =  .1165190029271822839255174533742327183599;
    A(9,8) = -.3401367187500000000000000000000000000000e-1;

    A(10,1) =  .4836757646340646986611287718844085773549e-1;
    A(10,2) =  0.;
    A(10,3) =  0.;
    A(10,4) =  0.;
    A(10,5) =  0.;
    A(10,6) =  .3928989925676163974333190042057047002852e-1;
    A(10,7) =  .1054740945890344608263649267140088017604;
    A(10,8) = -.2143865284648312665982642293830533996214e-1;
    A(10,9) = -.1041229174627194437759832813847147895623;

    A(11,1) = -.2664561487201478635337289243849737340534e-1;
    A(11,2) =  0.;
    A(11,3) =  0.;
    A(11,4) =  0.;
    A(11,5) =  0.;
    A(11,6) =  .3333333333333333333333333333333333333333e-1;
    A(11,7) = -.1631072244872467239162704487554706387141;
    A(11,8) =  .3396081684127761199487954930015522928244e-1;
    A(11,9) =  .1572319413814626097110769806810024118077;
    A(11,10) =  .2152267478031879552303534778794770376960;

    A(12,1) =  .3689009248708622334786359863227633989718e-1;
    A(12,2) =  0.;
    A(12,3) =  0.;
    A(12,4) =  0.;
    A(12,5) =  0.;
    A(12,6) = -.1465181576725542928653609891758501156785;
    A(12,7) =  .2242577768172024345345469822625833796001;
    A(12,8) =  .2294405717066072637090897902753790803034e-1;
    A(12,9) = -.35850052905728761357394424889330334334e-2;
    A(12,10) =  .8669223316444385506869203619044453906053e-1;
    A(12,11) =  .4383840651968337846196219974168630120572;

    A(13,1) = -.4866012215113340846662212357570395295088;
    A(13,2) =  0.;
    A(13,3) =  0.;
    A(13,4) =  0.;;
    A(13,5) =  0.;
    A(13,6) = -6.304602650282852990657772792012007122988;
    A(13,7) = -.281245618289472564778284183790118418111;
    A(13,8) = -2.679019236219849057687906597489223155566;
    A(13,9) =  .518815663924157511565311164615012522024;
    A(13,10) =  1.365353187603341710683633635235238678626;
    A(13,11) =  5.885091088503946585721274891680604830712;
    A(13,12) =  2.802808786272062889819965117517532194812;

    A(14,1) =  .4185367457753471441471025246471931649633;
    A(14,2) =  0.;
    A(14,3) =  0.;
    A(14,4) =  0.;
    A(14,5) =  0.;
    A(14,6) =  6.724547581906459363100870806514855026676;
    A(14,7) = -.425444280164611790606983409697113064616;
    A(14,8) =  3.343279153001265577811816947557982637749;
    A(14,9) =  .617081663117537759528421117507709784737;
    A(14,10) = -.929966123939932833937749523988800852013;
    A(14,11) = -6.099948804751010722472962837945508844846;
    A(14,12) = -3.002206187889399044804158084895173690015;
    A(14,13) =  .2553202529443445472336424602988558373637;

    A(15,1) = -.779374086122884664644623040843840506343;
    A(15,2) =  0.;
    A(15,3) =  0.;
    A(15,4) =  0.;
    A(15,5) =  0.;
    A(15,6) = -13.93734253810777678786523664804936051203;
    A(15,7) =  1.252048853379357320949735183924200895136;
    A(15,8) = -14.69150040801686878191527989293072091588;
    A(15,9) = -.494705058533141685655191992136962873577;
    A(15,10) =  2.242974909146236657906984549543692874755;
    A(15,11) =  13.36789380382864375813864978592679139881;
    A(15,12) =  14.39665048665068644512236935340272139005;
    A(15,13) = -.7975813331776800379127866056663258667437;
    A(15,14) =  .4409353709534277758753793068298041158235;

    A(16,1) =  2.058051337466886442151242368989994043993;
    A(16,2) =  0.;
    A(16,3) =  0.;
    A(16,4) =  0.;
    A(16,5) =  0.;
    A(16,6) =  22.35793772796803295519317565842520212899;
    A(16,7) =  .90949810997556332745009198137971890783;
    A(16,8) =  35.89110098240264104710550686568482456493;
    A(16,9) = -3.442515027624453437985000403608480262211;
    A(16,10) = -4.865481358036368826566013387928704014496;
    A(16,11) = -18.90980381354342625688427480879773032857;
    A(16,12) = -34.26354448030451782929251177395134170515;
    A(16,13) =  1.264756521695642578827783499806516664686;
    A(16,14) =  0.;
    A(16,15) =  0.;

    b(1) =  .1461197685842315252051541915018784713459e-1;
    b(2) =  0.;
    b(3) =  0.;
    b(4) =  0.;
    b(5) =  0.;
    b(6) =  0.;
    b(7) =  0.;
    b(8) = -.3915211862331339089410228267288242030810;
    b(9) =  .2310932500289506415909675644868993669908;
    b(10) =  .1274766769992852382560589467488989175618;
    b(11) =  .2246434176204157731566981937082069688984;
    b(12) =  .5684352689748512932705226972873692126743;
    b(13) =  .5825871557215827200814768021863420902155e-1;
    b(14) =  .1364317403482215641609022744494239843327;
    b(15) =  .3057013983082797397721005067920369646664e-1;
    b(16) =  0.;    
    
  otherwise
    error('invalid integrator option');
end