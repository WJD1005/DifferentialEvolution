

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outcome=jSO(func,totaltime,problem_size)

%     clc;
%     clear all;
%     format long;
%     format compact;
%     'SCSS_jSO'

%     func=15;
%     totaltime=8;
%     problem_size=30;
    
    rand('seed', sum(100 * clock));
    
    D = problem_size;
    pop_size = round(sqrt(D) * log(D) * 25);
    max_pop_size = pop_size;
    min_pop_size = 4;
    Max_FES = D * 10000;
    val_2_reach = 10^(-8);
    lu = [-100 * ones(1, D); 100 * ones(1, D)];

for func = func
  optimum = func * 100.0;

  %% Record the best results
  outcome = []; 

  fprintf('\n-------------------------------------------------------\n')
  fprintf('Function = %d, Dimension size = %d\n', func, D)
  
%   parfor run_id = 1 : 3
  parfor run_id = 1 : totaltime

    %% Initialize the main population
    popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, D) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
    pop = popold; % the old population becomes the current population

    fitness = cec17_func(pop',func);
    objfitness = fitness';

    popsize = pop_size;
    arc_rate = 1;
    arc_size = round(popsize * arc_rate);
    memory_size = 5;
    memory_pos = 1;
    memory_sf = 0.5 * ones(1, memory_size);
    memory_cr = 0.8 * ones(1, memory_size);

    pop_sf = zeros(1, popsize);
    pop_cr = zeros(1, popsize);
    g_p_best_rate = 0.25;
    p_best_rate = g_p_best_rate;
    p_num = round(popsize * p_best_rate);

    archive = [];
    arc_ind_count = 0;

    [bsf_fitness, minindex] = min(objfitness);
    bsf_solution = pop(minindex, :);
    FES = popsize;
    gen = 0;

    %% main loop
    while FES < Max_FES
    gen = gen + 1;
    child = zeros(popsize, D);
    [temp_fit, sorted_array] = sort(objfitness, 'ascend');

    for i = 1:popsize
        rmz = ceil(rand * memory_size);
        if rmz == memory_size
            mu_sf = 0.9; mu_cr = 0.9;
        else
            mu_sf = memory_sf(rmz);
            mu_cr = memory_cr(rmz);
        end
        if mu_cr < 0
            pop_cr(i) = 0;
        else
            pop_cr(i) = mu_cr + 0.1 * sqrt(-2 * log(rand)) * sin(2 * pi * rand); 
            if pop_cr(i) > 1
                pop_cr(i) = 1;
            elseif pop_cr(i) < 0
                pop_cr(i) = 0;
            end
        end
        if FES < 0.25 * Max_FES && pop_cr(i) < 0.7
            pop_cr(i) = 0.7;
        end
        if FES < 0.5 * Max_FES && pop_cr(i) < 0.6
            pop_cr(i) = 0.6;
        end
        
        pop_sf(i) = mu_sf + 0.1 * tan(pi * (rand - 0.5)); 
        while pop_sf(i) <= 0
            pop_sf(i) = mu_sf + 0.1 * tan(pi * (rand - 0.5));  
        end
        
        if pop_sf(i) > 1
            pop_sf(i) = 1;
        end
        if FES < 0.6 * Max_FES && pop_sf(i) > 0.7
            pop_sf(i) = 0.7;
        end
        
        p_best_ind = sorted_array(ceil(rand * p_num));
        while FES < 0.5 * Max_FES && p_best_ind == i
            p_best_ind = sorted_array(ceil(rand * p_num));   
        end

        jF = pop_sf(i);
        if FES < 0.2 * Max_FES
            jF = jF * 0.7;
        elseif FES < 0.4 * Max_FES
            jF = jF * 0.8;
        else
            jF = jF * 1.2;
        end
        r1 = ceil(rand * popsize);
        while r1 == i
            r1 = ceil(rand * popsize);
        end
        r2 = ceil(rand * (popsize + arc_ind_count));
        while r2 == i || r2 == r1
            r2 = ceil(rand * (popsize + arc_ind_count));
        end
        random_variable = ceil(rand * D);
        if r2 > popsize
            r2 = r2 - popsize;
            rrr = rand(1, D);
            from_el = rrr < pop_cr(i);
            from_el(random_variable) = 1;
            child(i, :) = from_el .* (pop(i, :) + jF * (pop(p_best_ind, :) - pop(i, :)) + pop_sf(i) * (pop(r1, :) - archive(r2, :))) + (1 - from_el) .* (pop(i, :)); 
        else
            rrr = rand(1, D);
            from_el = rrr < pop_cr(i);
            from_el(random_variable) = 1;
            child(i, :) = from_el .* (pop(i, :) + jF * (pop(p_best_ind, :) - pop(i, :)) + pop_sf(i) * (pop(r1, :) - pop(r2, :))) + (1 - from_el) .* (pop(i, :));
        end
        %---repair boundary
        less_min = child(i, :) < lu(1, :);  
        child(i, :) = (1 - less_min) .* child(i, :) + less_min .* ((pop(i, :) + lu(1, :)) / 2);
        more_max = child(i, :) > lu(2, :);
        child(i, :) = (1 - more_max) .* child(i, :) + more_max .* ((pop(i, :) + lu(2, :)) / 2);
        %-------------------  
    end  %end popsize ,end for
    

    
     child_fitness = zeros(1,popsize);
    for i = 1:popsize
      child_fitness(i) = cec17_func(child(i,:)',func);  
    end
    
    child_fitness = child_fitness';
    FES = FES + popsize;
    if FES < Max_FES
        [min_child_fitness, mindex] = min(child_fitness);
        if min_child_fitness < bsf_fitness
            bsf_fitness = min_child_fitness;
            bsf_solution = child(mindex(1), :);
        end
    else
        elp = Max_FES - (FES - popsize);
        elp_child_fitness = child_fitness(1:elp);
        elp_child = child(1:elp, :);
        [min_elp_child_fitness, mindex] = min(elp_child_fitness);
        if min_elp_child_fitness < bsf_fitness
            bsf_fitness = min_elp_child_fitness;
            bsf_solution = elp_child(mindex, :);
        end
        break;
    end
    %------------Add individuals to archive----------------
    equal_fit = (objfitness == child_fitness);
    pop(equal_fit, :) = child(equal_fit, :);
    
    less_fit = (child_fitness < objfitness);
    dif_fitness = abs(objfitness(less_fit) - child_fitness(less_fit));
    objfitness(less_fit) = child_fitness(less_fit);
    success_sf = pop_sf(less_fit);
    success_cr = pop_cr(less_fit);
   if arc_size > 1
    addinto_archive = pop(less_fit, :);  
    addinto_archive_num = size(addinto_archive, 1); 
    if (arc_ind_count + addinto_archive_num) <= arc_size
        archive = [archive; addinto_archive];
        arc_ind_count = arc_ind_count + addinto_archive_num;
    else
        releft_num = addinto_archive_num + arc_ind_count - arc_size;
        into_num = arc_size - arc_ind_count;
        archive = [archive; addinto_archive(1:into_num, :)];
        arc_ind_count = arc_size;
        addinto_archive(1:into_num, :) = [];
        for jj = 1:releft_num
            r11 = ceil(rand * arc_size);
            archive(r11, :) = addinto_archive(jj, :);
        end
    end
   end %end if arc_size > 1
    pop(less_fit, :) = child(less_fit, :);
    %-----------------------------------------------
    %old_num_success_params = num_success_params;
    num_success_params = length(success_sf);
    if num_success_params > 0
        old_sf = memory_sf(memory_pos);
        old_cr = memory_cr(memory_pos);
        memory_sf(memory_pos) = 0;
        memory_cr(memory_pos) = 0;

        sumA = sum(dif_fitness);
        weightA = dif_fitness / sumA;  %col vector
        all_sf_product = weightA' .* success_sf .* success_sf; %succexx_sf is row vector
        tem_sf_product = weightA' .* success_sf;
        memory_sf(memory_pos) = sum(all_sf_product);
        temp_sum_sf = sum(tem_sf_product);
        all_cr_product = weightA' .* success_cr .* success_cr;
        tem_cr_product = weightA' .* success_cr;
        memory_cr(memory_pos) = sum(all_cr_product);
        temp_sum_cr = sum(tem_cr_product);
        
        memory_sf(memory_pos) = memory_sf(memory_pos) / temp_sum_sf;
        if temp_sum_cr == 0 || memory_cr(memory_pos) == -1
            memory_cr(memory_pos) = -1;
        else
            memory_cr(memory_pos) = memory_cr(memory_pos) / temp_sum_cr;
        end
        memory_sf(memory_pos) = (memory_sf(memory_pos) + old_sf) / 2;
        memory_cr(memory_pos) = (memory_cr(memory_pos) + old_cr) / 2;
        
        memory_pos = memory_pos + 1;
        if memory_pos > memory_size
            memory_pos = 1;
        end
    end
    %-----------Calculate the next population size-----------------------
    plan_pop_size = round((min_pop_size - max_pop_size) / (Max_FES) * FES + max_pop_size);
    if popsize > plan_pop_size
        reduction_ind_num = popsize - plan_pop_size;
        if (popsize - reduction_ind_num) < min_pop_size
            reduction_ind_num = popsize - min_pop_size;
        end
        [objfitness, sindexx] = sort(objfitness, 'ascend');
        objfitness((popsize - reduction_ind_num + 1) : popsize) = [];
        pop = pop(sindexx, :);
        pop((popsize - reduction_ind_num + 1) : popsize, :) = [];
        popsize = popsize - reduction_ind_num;
        arc_size = popsize * arc_rate;
        if arc_ind_count > arc_size
            arc_ind_count = arc_size;
        end
        p_best_rate = g_p_best_rate * (1 - 0.5 * FES / Max_FES);
        p_num = round(popsize * p_best_rate);
        if p_num <=1
            p_num = 2;
        end 
    end  
    end

    bsf_error_val = min(objfitness) - optimum;
    if bsf_error_val < val_2_reach
        bsf_error_val = 0;
     end

%     fprintf('jSO_%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
    outcome = [outcome bsf_error_val];   
     
  end %% end 1 run
  
  fprintf('\n')
  fprintf('jSO_mean error value = %1.8e, std = %1.8e\n', mean(outcome), std(outcome))
% %    if problem_size == 30
% %   rowid = sprintf('A%d', func);
% %   xlswrite('jSO.xls',sort(outcome),'D30',rowid); 
% %   xlswrite('jSO.xls',[mean(outcome),std(outcome)],'D30mean&std',rowid); 
% %   elseif problem_size == 50
% %   rowid = sprintf('A%d', func);
% %   xlswrite('jSO.xls',sort(outcome),'D50',rowid); 
% %   xlswrite('jSO.xls',[mean(outcome),std(outcome)],'D50mean&std',rowid);
% %    elseif problem_size == 100
% %   rowid = sprintf('A%d', func);
% %   xlswrite('jSO.xls',sort(outcome),'D100',rowid); 
% %   xlswrite('jSO.xls',[mean(outcome),std(outcome)],'D100mean&std',rowid);
% %   end  
end %% end 1 function run
end



