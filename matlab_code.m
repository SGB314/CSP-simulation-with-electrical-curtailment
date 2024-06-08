clear all, clc, close all

tic
%% Questions

year = input('Introduce un año (2016-2021): ');
curtailments = input('¿Quieres cálculo con curtailments (si/no)?: ', 's');
graphics = input('¿Quieres generar las gráficas diarias (si/no)?: ', 's');

%% Load data


horas = 1:(24*365); % 24h/dia x 365 días = 8760h
diezminutales = linspace(0,8760,52561); % 6 diezminutos/hora x 24h x 365 = 52560
step = 6; % 6 diezminutos en 1 hora

%Irradiation data + interpolation (conversion 1h => 10min)

DataSAM_1h = load(strcat(pwd,'\input\','\sam_results\','results.txt')); %SAM: Thermal energy incident on receiver (MWh)
DataSAM_10min = interp1(horas,DataSAM_1h,diezminutales,"pchip")';
time_stamp = load(strcat(pwd,'\input\','\sam_results\','time_stamp.txt'));

contador = 0;


%Surplus data + interpolation (conversion 1h => 10min)


% For Leap Years it eliminate last row.

if (year == 2016 || year == 2020)

    wind_1h = load(fullfile(pwd, 'input', 'surplus_electricity', sprintf('%u_wind_1h.txt', year))); % MWh
    wind_1h_wo_onerow = wind_1h(1:end-1, :);
    wind_10min = abs(interp1(horas, reshape(wind_1h_wo_onerow.', [], 1), diezminutales, 'pchip'))';
    null_wind = find(isnan(wind_10min));
    wind_10min(null_wind,1)=0;

    fot_1h = load(fullfile(pwd, 'input', 'surplus_electricity', sprintf('%u_fot_1h.txt', year))); % MWh
    fot_1h_wo_onerow = fot_1h(1:end-1, :);
    fot_10min = abs(interp1(horas, reshape(fot_1h_wo_onerow.', [], 1), diezminutales, 'pchip'))';
    null_fot = find(isnan(fot_10min));
    fot_10min(null_fot,1)=0;

    therm_1h = load(fullfile(pwd, 'input', 'surplus_electricity', sprintf('%u_therm_1h.txt', year))); % MWh
    therm_1h_wo_onerow = therm_1h(1:end-1, :);
    therm_10min = abs(interp1(horas, reshape(therm_1h_wo_onerow.', [], 1), diezminutales, 'pchip'))';
    null_therm = find(isnan(therm_10min));
    therm_10min(null_therm,1)=0;

elseif (year == 2017 || year == 2018 || year == 2019 || year == 2021)
    wind_1h = load(fullfile(pwd, 'input', 'surplus_electricity', sprintf('%u_wind_1h.txt', year))); % MWh
    wind_10min = abs(interp1(horas,reshape(wind_1h.',[],1),diezminutales,"pchip"))';
    null_wind = find(isnan(wind_10min));
    wind_10min(null_wind,1)=0;

    fot_1h = load(fullfile(pwd, 'input', 'surplus_electricity', sprintf('%u_fot_1h.txt', year))); % MWh
    fot_10min = abs(interp1(horas,reshape(fot_1h.',[],1),diezminutales,"pchip"))';
    null_fot = find(isnan(fot_10min));
    fot_10min(null_fot,1)=0;

    therm_1h = load(fullfile(pwd, 'input', 'surplus_electricity', sprintf('%u_therm_1h.txt', year))); % MWh
    therm_10min = abs(interp1(horas,reshape(therm_1h.',[],1),diezminutales,"pchip"))';
    null_therm = find(isnan(therm_10min));
    therm_10min(null_therm,1)=0;
else
    year = input('Introduce un año válido (2016-2021): ');
end



%% GENERAL PARAMETERS

% SCENARIO I/II


% eff_rec = 0.885;         % Receiver efficiency
% eff_hx = 0.85;       % Heat Exchanger Efficiency
% eff_pb = 0.412;          % Power Block Efficiency
% Wnet_nom = 19.9; %MW     % Nominal power of the power block
% storage_losses = (pi*((0.5*14.8774)^2)*72.70 + pi*((0.5*14.8774)^2)*61 + 2*pi*(0.5*14.8774)*20*79.13)*10^-6; %MW ref: https://doi.org/10.1016/j.solener.2016.06.030
% eff_transformer = 0.975;
% eff_resistance = 0.975;
% eff_carnot = eff_transformer * eff_resistance; % Scenario II


% SCENARIO III/IV


eff_rec = 0.795;         % Receiver efficiency
eff_hx = 0.96;       % Heat Exchanger Efficiency
eff_pb = 0.454;          % Power Block Efficiency
Wnet_nom = 19.9; %MW     % Nominal power of the power block
storage_losses = (pi*((0.5*14.8774)^2)*72.70 + pi*((0.5*14.8774)^2)*61 + 2*pi*(0.5*14.8774)*20*79.13)*10^-6; %MW ref: https://doi.org/10.1016/j.solener.2016.06.030
eff_transformer = 0.975;
eff_resistance = 0.975;
eff_carnot = eff_transformer * eff_resistance * eff_hx; % Scenario IV



%% INITIAL DATA

neg_ind = find(isnan(DataSAM_10min));
DataSAM_10min(neg_ind,1)=0;
Wnet(1,1) = 0;
used(1,1) = 0;
tot_vector = [];
Sum_surplus = 0;
Sum_storage_cold = 0;
Sum_We = 0;
partial = 0;


%% STORAGE DATA

storage_max = 15 * (Wnet_nom/(eff_hx*eff_pb)) + (5/100) * 15 * (Wnet_nom/(eff_hx*eff_pb)); %MWh
storage_min = (5/100) * storage_max; %MWh
storage_produce = 4 * storage_min;
storage_ini = storage_min;



resistance = 15; % MW



t_start_shutdown = 1; %Must be a multiple of "step"
variation = 1/(t_start_shutdown * step);



%% Calculation of dinamic power block efficiency


load_PB = [0 0.60 0.85 1];
efficiency_PB = (eff_pb/0.37698) * [0 0.35958 0.37167 0.37698]; %REF:https://doi.org/10.18186/thermal.298611

degree = 3;

coefficients = polyfit(load_PB,efficiency_PB,degree);


%% START LOOP


for n_Dia = 1:365

    ini_Dia_TRNSYS = (n_Dia - 1) * 24 * step + 1; % cada día abarca 144 valores porque son diezminutales
    fin_Dia_TRNSYS = n_Dia * 24 * step;


    storage(ini_Dia_TRNSYS,1) = storage_ini;




    for k = ini_Dia_TRNSYS + 1 : fin_Dia_TRNSYS + 1




        if strcmpi(curtailments, 'si') || strcmpi(curtailments, 'sí')



            sun(k,1) = (((DataSAM_10min(k,1) + DataSAM_10min(k-1,1)) * eff_rec / 2 *(1/step)));
            surplus_electricity(k,1) = (eff_carnot/(2*step)) * ((wind_10min(k) + wind_10min(k-1)) + (fot_10min(k) + fot_10min(k-1)) + (therm_10min(k) + therm_10min(k-1)));

            % STORAGE

            if surplus_electricity(k,1) >= resistance
                surplus_electricity(k,1) = resistance;
            else
                surplus_electricity(k,1) = (eff_carnot/(2*step)) * ((wind_10min(k) + wind_10min(k-1)) + (fot_10min(k) + fot_10min(k-1)) + (therm_10min(k) + therm_10min(k-1)));
            end

            storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses + surplus_electricity(k,1);



            if (storage(k,1) > storage_min) && (storage(k,1) < storage_max)
                if (storage_max - storage (k,1)) >= surplus_electricity(k,1)
                    if surplus_electricity(k,1) >= resistance
                        surplus_electricity(k,1) = resistance;
                        surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                        storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses + surplus_electricity(k,1);
                    elseif surplus_electricity(k,1) == 0
                        storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses;
                        surplus_electricity(k,1) = 0;
                        surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                    else
                        surplus_electricity(k,1) = (eff_carnot/(2*step)) * ((wind_10min(k) + wind_10min(k-1)) + (fot_10min(k) + fot_10min(k-1)) + (therm_10min(k) + therm_10min(k-1)));
                        surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                        storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses + surplus_electricity(k,1);
                    end

                else
                    if surplus_electricity(k,1) > 0
                        surplus_electricity(k,1) = (storage_max - storage (k,1));
                        surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                        storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses + surplus_electricity(k,1);
                    else
                        storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses;
                        surplus_electricity(k,1) = 0;
                        surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                    end
                end
            elseif storage(k,1) >= storage_max
                storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses;
                if storage(k,1) >= storage_max
                    storage(k,1) = storage_max;
                    surplus_electricity(k,1) = 0;
                    surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                else
                    if (storage_max - storage (k,1)) > surplus_electricity(k,1)
                        if surplus_electricity(k,1) >= resistance
                            surplus_electricity(k,1) = resistance;
                            surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                            storage(k,1) = storage(k-1,1) - storage_losses + surplus_electricity(k,1) + sun(k,1);
                        elseif surplus_electricity(k,1) == 0
                            storage(k,1) = storage(k-1,1) - storage_losses + sun(k,1);
                            surplus_electricity(k,1) = 0;
                            surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                        else
                            surplus_electricity(k,1) = (eff_carnot/(2*step)) * ((wind_10min(k) + wind_10min(k-1)) + (fot_10min(k) + fot_10min(k-1)) + (therm_10min(k) + therm_10min(k-1)));
                            surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                            storage(k,1) = storage(k-1,1) - storage_losses + surplus_electricity(k,1) + sun(k,1);
                        end

                    else
                        if surplus_electricity(k,1) > 0
                            surplus_electricity(k,1) = (storage_max - storage (k,1));
                            surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                            storage(k,1) = storage(k-1,1) - storage_losses + surplus_electricity(k,1) + sun(k,1);
                        else
                            storage(k,1) = storage(k-1,1) - storage_losses + sun(k,1);
                            surplus_electricity(k,1) = 0;
                            surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                        end
                    end
                end

            else
                contador = contador + 1;
                if surplus_electricity(k,1) >= resistance
                    surplus_electricity(k,1) = resistance;
                    surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                    storage(k,1) = storage(k-1,1) + sun(k,1) + surplus_electricity(k,1) - storage_losses;
                elseif surplus_electricity(k,1) == 0
                    storage(k,1) = storage(k-1,1) + sun(k,1);
                    surplus_electricity(k,1) = 0;
                    surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                else
                    surplus_electricity(k,1) = (eff_carnot/(2*step)) * ((wind_10min(k) + wind_10min(k-1)) + (fot_10min(k) + fot_10min(k-1)) + (therm_10min(k) + therm_10min(k-1)));
                    surplus_electricity_gross(k,1) = surplus_electricity(k,1) / eff_carnot;
                    storage(k,1) = storage(k-1,1) + sun(k,1) + surplus_electricity(k,1) - storage_losses;
                end

            end



            % PRODUCTION

            if partial <= 0.01 % PB = OFF
                if storage(k,1) <= storage_produce + storage_min
                    partial = 0;
                    Wnet(k,1) = 0;
                    used(k,1) = 0;
                    storage(k,1) = storage(k,1) - used(k,1);

                else
                    partial = partial + variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);

                end
            elseif partial >= 1 - 0.01 %PB = 100%
                if storage(k,1) >= (storage_min + 0.5*storage_produce)
                    partial = 1;
                    Wnet(k,1) = Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);

                else
                    partial = partial - variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);

                end
            else %during start/shutdown
                if storage(k,1) >= (storage_min + storage_produce - 0.5*storage(k-1,1))
                    partial = partial + variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);


                else
                    partial = partial - variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);

                end
            end


        elseif strcmpi(curtailments, 'no')



            sun(k,1) = (((DataSAM_10min(k,1) + DataSAM_10min(k-1,1)) * eff_rec / 2 *(1/step)));

            storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses;
            surplus_electricity(k,1) = 0;
            surplus_electricity_gross(k,1) = 0;

            % STORAGE

            if (storage(k,1) > storage_min) && (storage(k,1) < storage_max)
                storage(k,1) = storage(k-1,1) + sun(k,1) - storage_losses;
            elseif storage(k,1) >= storage_max
                storage(k,1) = storage_max;
            else
                storage(k,1) = storage(k-1,1) + sun(k,1);
                contador = contador + 1;
            end



            % PRODUCTION


            if partial <= 0.01 % PB=OFF
                if storage(k,1) <= storage_produce + storage_min
                    partial = 0;
                    Wnet(k,1) = 0;
                    used(k,1) = 0;
                    storage(k,1) = storage(k,1) - used(k,1);
                    storage_cold(k,1) = storage_max - storage(k,1);

                else
                    partial = partial + variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);
                    storage_cold(k,1) = storage_max - storage(k,1);

                end
            elseif partial >= 1 - 0.01 %PB = 100%
                if storage(k,1) >= (storage_min + 0.5*storage_produce)
                    partial = 1;
                    Wnet(k,1) = Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);
                    storage_cold(k,1) = storage_max - storage(k,1);

                else
                    partial = partial - variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);
                    storage_cold(k,1) = storage_max - storage(k,1);

                end
            else %during start/shutdown
                if storage(k,1) >= (storage_min + storage_produce)
                    partial = partial + variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);
                    storage_cold(k,1) = storage_max - storage(k,1);


                else
                    partial = partial - variation;
                    Wnet(k,1) = partial * Wnet_nom;
                    used(k,1) = (Wnet(k,1)/(polyval(coefficients,(Wnet(k,1)/Wnet_nom)) * eff_hx))/step;
                    storage(k,1) = storage(k,1) - used(k,1);
                    storage_cold(k,1) = storage_max - storage(k,1);

                end
            end
        else
            curtailments = input('Respuesta no válida. Por favor, responde con "si" o "no":', 's');
        end
    end

    storage_ini = storage(k,1);


    for j = ini_Dia_TRNSYS : fin_Dia_TRNSYS
        Sum_We = Sum_We + (Wnet (j) / step);
    end

    datos_figuras{n_Dia} = {time_stamp(ini_Dia_TRNSYS:fin_Dia_TRNSYS, 1), DataSAM_10min(ini_Dia_TRNSYS:fin_Dia_TRNSYS, 1), Wnet(ini_Dia_TRNSYS:fin_Dia_TRNSYS, 1), storage(ini_Dia_TRNSYS:fin_Dia_TRNSYS, 1), surplus_electricity(ini_Dia_TRNSYS:fin_Dia_TRNSYS, 1)};

end

fprintf('The total electricity production is %d MWh\n', Sum_We);



%%FIGURES

if strcmpi(graphics, 'si') && strcmpi(curtailments, 'si')

    for n_Dia = 1:365



        % Extrae los datos relevantes para esta figura
        datos = datos_figuras{n_Dia};
        x_values = datos{1};
        DataSAM_10min_values = datos{2};
        Wnet_values = datos{3};
        storage_values = datos{4};
        surplus_values = datos{5};





        % Dibuja la figura
        figure('Color', [1 1 1], 'Position', [2000 100 600 400], 'Name', strcat('Dia', num2str(n_Dia)));
        hold on; box on; grid on;

        date = datetime(2021,1,n_Dia);
        date_str = [num2str(day(date)), month(date, 'long')];
        title([date_str])





        xlim([x_values(1), x_values(end)]);
        xlabel('Horas');

        % Trazar la radiación incidente en el eje izquierdo

        yyaxis left;
        ylabel('Potencia [MW]');
        ylim([0 175]);  % Ajusta el límite superior según los datos
        yticks(0:25:175);
        plot(x_values, DataSAM_10min_values, 'LineStyle', '-', 'Marker', '.', 'Color', 'magenta', 'DisplayName', 'Irradiation');


        % Trazar los datos de Wnet en el eje izquierdo

        plot(x_values, Wnet_values,'LineStyle','-','Color','green','Marker','.', 'DisplayName','Power block');
        plot(x_values, surplus_values,"LineWidth",1,"Color",'[0.9290 0.6940 0.1250]','LineStyle','-','Marker','.');


        % Trazar los datos de storage en el eje derecho
        yyaxis right;
        ylabel('Almacenamiento [MWh]');
        ylim([0 950]);
        yticks(0:100:950);
        plot(x_values, storage_values, 'LineStyle', '-', 'Marker', '.', 'DisplayName', 'Storage');
        line(xlim, [storage_min storage_min],'Color','k','LineStyle', '--');
        line(xlim, [storage_max storage_max],'Color','k', 'LineStyle', '--');
        legend('Irradiación','BP','Vertidos','Almacenamiento','Almacenamiento mín-máx','Location','northoutside','Orientation','vertical')
        grid on



        xticks(1:1:52560);
        ticks = { '', '', '3h', '', '', '6h', '', '', '9h', '', '', '12h', '', '', '15h', '', '', '18h', '', '', '21h', '', '', ''};
        hours = repmat(ticks, 1, n_Dia);
        xticklabels(hours);



        % Guardar la figura
        saveas(gcf, strcat(pwd, '\figures\with surplus energy\', num2str(year), '\', num2str(year), '_Day_', num2str(n_Dia), '.png'));
        close(gcf);
    end
elseif strcmpi(graphics, 'si') && strcmpi(curtailments, 'no')
    for n_Dia = 1:365



        % Extrae los datos relevantes para esta figura
        datos = datos_figuras{n_Dia};
        x_values = datos{1};
        DataSAM_10min_values = datos{2};
        Wnet_values = datos{3};
        storage_values = datos{4};





        % Dibuja la figura
        figure('Color', [1 1 1], 'Position', [2000 100 600 400], 'Name', strcat('Dia', num2str(n_Dia)));
        hold on; box on; grid on;

        date = datetime(2021,1,n_Dia);
        date_str = [num2str(day(date)), month(date, 'long')];
        title([date_str])


        xlim([x_values(1), x_values(end)]);
        xlabel('Horas');

        % Trazar la radiación incidente en el eje izquierdo

        yyaxis left;
        ylabel('Potencia [MW]');
        ylim([0 175]);  % Ajusta el límite superior según los datos
        yticks(0:25:175);
        plot(x_values, DataSAM_10min_values, 'LineStyle', '-', 'Marker', '.', 'DisplayName', 'Radiation','Color','magenta');


        % Trazar los datos de Wnet en el eje izquierdo

        plot(x_values, Wnet_values,'LineStyle','-','Marker','.', 'DisplayName','Power block','Color','green');


        % Trazar los datos de storage en el eje derecho
        yyaxis right;
        ylabel('Almacenamiento [MWh]');
        ylim([0 950]);
        plot(x_values, storage_values, 'LineStyle', '-', 'Marker', '.', 'DisplayName', 'Storage');

        line(xlim, [storage_min storage_min],'Color','k','LineStyle', '--');
        line(xlim, [storage_max storage_max],'Color','k', 'LineStyle', '--');
        legend('Irradiación','BP','Almacenamiento','Almacenamiento mín-máx','Location','northoutside','Orientation','vertical')



        xticks(1:1:52560);
        ticks = { '', '', '3h', '', '', '6h', '', '', '9h', '', '', '12h', '', '', '15h', '', '', '18h', '', '', '21h', '', '', ''};
        hours = repmat(ticks, 1, n_Dia);
        xticklabels(hours);



        % Guardar la figura
        saveas(gcf, strcat(pwd, '\figures\without surplus energy\', 'Day_', num2str(n_Dia), '.png'));
        close(gcf);
    end
end


toc
