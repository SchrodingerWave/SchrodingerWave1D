classdef SchrodingerWave < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        ExportPotentialMenu_3           matlab.ui.container.Menu
        ExportEigenvaluescsvMenu_2      matlab.ui.container.Menu
        ExportEigenvectorscsvMenu       matlab.ui.container.Menu
        SettingMenu                     matlab.ui.container.Menu
        ThemeMenu                       matlab.ui.container.Menu
        WhiteMenu                       matlab.ui.container.Menu
        BlueMenu                        matlab.ui.container.Menu
        DarkMenu                        matlab.ui.container.Menu
        DefaultMenu                     matlab.ui.container.Menu
        OwlsMenu                        matlab.ui.container.Menu
        ShowmemoreMenu                  matlab.ui.container.Menu
        FlythemoffMenu                  matlab.ui.container.Menu
        ResolutionMenu                  matlab.ui.container.Menu
        SettoHighMenu                   matlab.ui.container.Menu
        SettoMediumMenu                 matlab.ui.container.Menu
        SettoLowMenu                    matlab.ui.container.Menu
        InfoonResolutionMenu            matlab.ui.container.Menu
        HiddenPotentialSettingsMenu     matlab.ui.container.Menu
        CoulmbicMenu                    matlab.ui.container.Menu
        SettingMenu_4                   matlab.ui.container.Menu
        HelpMenu_4                      matlab.ui.container.Menu
        HarmonicOscilatorMenu           matlab.ui.container.Menu
        SettingMenu_3                   matlab.ui.container.Menu
        HelpMenu_3                      matlab.ui.container.Menu
        MorseOscillatorMenu             matlab.ui.container.Menu
        SettingMenu_2                   matlab.ui.container.Menu
        HelpMenu_2                      matlab.ui.container.Menu
        AllowNegativeValuesofPotentialEnergyMenu  matlab.ui.container.Menu
        ShiftWavefunctionbyEnergyMenu   matlab.ui.container.Menu
        EnableorDisableMenu             matlab.ui.container.Menu
        ShowEnergySpectrumMenu          matlab.ui.container.Menu
        SelectandPlotEigenvaluesMenu    matlab.ui.container.Menu
        EnergySpectrumPlottingOptionsMenu  matlab.ui.container.Menu
        ChangePlotTypeMenu              matlab.ui.container.Menu
        NormalizationConstantMenu       matlab.ui.container.Menu
        ChangeNormalizationConstantMenu  matlab.ui.container.Menu
        HelpMenu                        matlab.ui.container.Menu
        DocumentationMenu               matlab.ui.container.Menu
        AboutMenu                       matlab.ui.container.Menu
        UIAxes                          matlab.ui.control.UIAxes
        RunButton                       matlab.ui.control.Button
        PeriodicBoundaryConditionDropDownLabel  matlab.ui.control.Label
        PeriodicBoundaryConditionDropDown  matlab.ui.control.DropDown
        PotentialEnergyPanel            matlab.ui.container.Panel
        FiniteSquareWellCheckBox        matlab.ui.control.CheckBox
        FiniteSquareBarrierCheckBox     matlab.ui.control.CheckBox
        CoulombPotentialCheckBox        matlab.ui.control.CheckBox
        StepPotentialCheckBox           matlab.ui.control.CheckBox
        HarmonicOscillatorCheckBox      matlab.ui.control.CheckBox
        MorseOscillatorCheckBox         matlab.ui.control.CheckBox
        RandomPotentialCheckBox         matlab.ui.control.CheckBox
        DoubleWellCheckBox              matlab.ui.control.CheckBox
        ImportedPotentialCheckBox       matlab.ui.control.CheckBox
        EigenstateVisualizationPanel    matlab.ui.container.Panel
        TabGroup                        matlab.ui.container.TabGroup
        SingleTab                       matlab.ui.container.Tab
        EigenStateQuantumNumberLabel    matlab.ui.control.Label
        EigenStateQuantumNumberSpinner  matlab.ui.control.Spinner
        Wavefunction2Label              matlab.ui.control.Label
        Wavefunction2Switch             matlab.ui.control.Switch
        MultiStateVisualizationmodeLabel  matlab.ui.control.Label
        MultiTab                        matlab.ui.container.Tab
        MultiSwitch                     matlab.ui.control.Switch
        StateSpinnerLabel               matlab.ui.control.Label
        StateSpinner                    matlab.ui.control.Spinner
        StateSpinner_2Label             matlab.ui.control.Label
        StateSpinner_2                  matlab.ui.control.Spinner
        StateSpinner_3Label             matlab.ui.control.Label
        StateSpinner_3                  matlab.ui.control.Spinner
        StateSpinner_4Label             matlab.ui.control.Label
        StateSpinner_4                  matlab.ui.control.Spinner
        Image                           matlab.ui.control.Image
        MassDropDownLabel               matlab.ui.control.Label
        MassDropDown                    matlab.ui.control.DropDown
        PotentialModificationPanel      matlab.ui.container.Panel
        TabGroup2                       matlab.ui.container.TabGroup
        WellTab                         matlab.ui.container.Tab
        WellDepthEditFieldLabel         matlab.ui.control.Label
        WellDepthEditField              matlab.ui.control.NumericEditField
        LeftWallPositionEditFieldLabel  matlab.ui.control.Label
        LeftWallPositionEditField       matlab.ui.control.NumericEditField
        RightWallPositionEditFieldLabel  matlab.ui.control.Label
        RightWallPositionEditField      matlab.ui.control.NumericEditField
        ResettoDefaultButton            matlab.ui.control.Button
        BarrierTab                      matlab.ui.container.Tab
        ResettoDefaultButton_2          matlab.ui.control.Button
        LeftWallPositionEditField_2Label  matlab.ui.control.Label
        LeftWallPositionEditField_2     matlab.ui.control.NumericEditField
        RightWallPositionEditField_2Label  matlab.ui.control.Label
        RightWallPositionEditField_2    matlab.ui.control.NumericEditField
        BarrierHeightEditFieldLabel     matlab.ui.control.Label
        BarrierHeightEditField          matlab.ui.control.NumericEditField
        StepTab                         matlab.ui.container.Tab
        WallPositionEditFieldLabel      matlab.ui.control.Label
        WallPositionEditField           matlab.ui.control.NumericEditField
        RightPotentialEditFieldLabel    matlab.ui.control.Label
        RightPotentialEditField         matlab.ui.control.NumericEditField
        LeftPotentialEditFieldLabel     matlab.ui.control.Label
        LeftPotentialEditField          matlab.ui.control.NumericEditField
        ResettoDefaultButton_3          matlab.ui.control.Button
        TabGroup3                       matlab.ui.container.TabGroup
        HarmonicTab                     matlab.ui.container.Tab
        ResettoDefaultButton_4          matlab.ui.control.Button
        EquPositionEditFieldLabel       matlab.ui.control.Label
        EquPositionEditField            matlab.ui.control.NumericEditField
        ForceConstantEditFieldLabel     matlab.ui.control.Label
        ForceConstantEditField          matlab.ui.control.NumericEditField
        MorseTab                        matlab.ui.container.Tab
        ResettoDefaultButton_5          matlab.ui.control.Button
        ForceConstantEditField_2Label   matlab.ui.control.Label
        ForceConstantEditField_2        matlab.ui.control.NumericEditField
        EquPositionEditField_2Label     matlab.ui.control.Label
        EquPositionEditField_2          matlab.ui.control.NumericEditField
        DissociationEnergyEditFieldLabel  matlab.ui.control.Label
        DissociationEnergyEditField     matlab.ui.control.NumericEditField
        RandomTab_2                     matlab.ui.container.Tab
        NoiseWeightEditFieldLabel       matlab.ui.control.Label
        NoiseWeightEditField            matlab.ui.control.NumericEditField
        ResettoDefaultButton_7          matlab.ui.control.Button
        TabGroup4                       matlab.ui.container.TabGroup
        DoubelWellTab                   matlab.ui.container.Tab
        PEcaxx1xx2xx3xx4Label           matlab.ui.control.Label
        aEditFieldLabel                 matlab.ui.control.Label
        aEditField                      matlab.ui.control.NumericEditField
        x1EditFieldLabel                matlab.ui.control.Label
        x1EditField                     matlab.ui.control.NumericEditField
        x2EditFieldLabel                matlab.ui.control.Label
        x2EditField                     matlab.ui.control.NumericEditField
        x3EditFieldLabel                matlab.ui.control.Label
        x3EditField                     matlab.ui.control.NumericEditField
        x4EditFieldLabel                matlab.ui.control.Label
        x4EditField                     matlab.ui.control.NumericEditField
        ResettoDefaultButton_6          matlab.ui.control.Button
        cEditFieldLabel                 matlab.ui.control.Label
        cEditField                      matlab.ui.control.NumericEditField
        ImportTab                       matlab.ui.container.Tab
        ImportPotentialButton           matlab.ui.control.StateButton
        IMPRTPETextArea                 matlab.ui.control.TextArea
        IMPRT_PELamp                    matlab.ui.control.Lamp
        UIAxes2                         matlab.ui.control.UIAxes
        Button                          matlab.ui.control.StateButton
        Image2                          matlab.ui.control.Image
        Image3                          matlab.ui.control.Image
        Image4                          matlab.ui.control.Image
    end

    
    properties (Access = private)
        setMinmum2Zero2Well=1; % Description
        IMP_TRUE    =0;
        PRVpath     ='';
        IPF         = [];
        IMPPE_File  ='';
        IMPPE_Path  ='';
        ExportFileName2='';ExportPathName2='';
        Finalpotential=[];
        Npoints=1000;
        StepSize=0.05;
        nNnTimes_morse=9;
        MaxLevelMorse=1.5;
        AddEigv2EigfPlot=1;
        EigVal=[];
        EigVec=[];
        CloumbicasymptotyicValue=2;
        NormalizationSwitch=1;
        NormalizingConstant=1;
        FirstParticleCharge=+1;
        SecondParticleCharge=-1;
        MaximumAllowedCloumbic_Harmonic=10; %Hartrees
        AllowNegatives=1;
        PlotType=1;
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.StateSpinner.BackgroundColor='k';
            app.StateSpinner_2.BackgroundColor='r';
            app.StateSpinner_3.BackgroundColor='b';
            app.StateSpinner_4.BackgroundColor='g';
            app.LeftWallPositionEditField.Enable='off';
            app.WellDepthEditField.Enable='off';
            app.RightWallPositionEditField.Enable='off';
            app.ResettoDefaultButton.Enable='off';
            app.LeftWallPositionEditField_2.Enable='off';
            app.RightWallPositionEditField_2.Enable='off';
            app.BarrierHeightEditField.Enable='off';
            app.ResettoDefaultButton_2.Enable='off';
            app.WallPositionEditField.Enable='off';
            app.RightPotentialEditField.Enable='off';
            app.LeftPotentialEditField.Enable='off';
            app.ResettoDefaultButton_3.Enable='off';
            app.EquPositionEditField.Enable='off';
            app.ForceConstantEditField.Enable='off';
            app.ResettoDefaultButton_4.Enable='off';
            app.ResettoDefaultButton_5.Enable='off';
            app.ForceConstantEditField_2.Enable='off';
            app.EquPositionEditField_2.Enable='off';
            app.DissociationEnergyEditField.Enable='off';
            app.ResettoDefaultButton_6.Enable='off';
            app.aEditField.Enable='off';
            app.x1EditField.Enable='off';
            app.x2EditField.Enable='off';
            app.x3EditField.Enable='off';
            app.x4EditField.Enable='off';
            app.PEcaxx1xx2xx3xx4Label.Enable='off';
            app.Image4.Visible='off';
            app.Image3.Visible='off';
            app.Image2.Visible='off';
            app.ResettoDefaultButton_7.Enable='off';
            app.NoiseWeightEditField.Enable='off';
            app.ImportPotentialButton.Enable='off';
            app.IMPRTPETextArea.Enable='off';
            app.IMPRT_PELamp.Enable='off';
            app.IMPRT_PELamp.Color='b';
            app.UIAxes2.Visible='off';
            app.PRVpath     =pwd;
        end

        % Value changed function: FiniteSquareWellCheckBox
        function FiniteSquareWellCheckBoxValueChanged(app, event)
            if app.FiniteSquareWellCheckBox.Value
                app.LeftWallPositionEditField.Enable='on';
                app.WellDepthEditField.Enable='on';
                app.RightWallPositionEditField.Enable='on';
                app.ResettoDefaultButton.Enable='on';
            else
                app.LeftWallPositionEditField.Enable='off';
                app.WellDepthEditField.Enable='off';
                app.RightWallPositionEditField.Enable='off';
                app.ResettoDefaultButton.Enable='off';
            end
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            cla(app.UIAxes)
            %% Define Constants

            hbar = 1;
            step_spacing =app.StepSize;
            DataPointsNumber = app.Npoints;
            if app.MassDropDown.Value=="Electron"
                mass=1;
            elseif app.MassDropDown.Value=="Proton"
                mass=1822;
            elseif app.MassDropDown.Value=="Helium"
                mass=7296;
            elseif app.MassDropDown.Value=="H2"
                mass=3660;
            end
            %%  Simulation Parameters

            if app.PeriodicBoundaryConditionDropDown.Value=="No"
                PerBoundCond = false;
            elseif app.PeriodicBoundaryConditionDropDown.Value=="Yes"
                PerBoundCond = true;
            end
            %%  Potential Energy
            FiniteSquareWell = false;
            SquareBarrier = false;
            CoulmbicInteractionPotential = false;
            StepFuncPotential = false;
            HarmonicOscillator = false;
            RandomNoise = false;
            morse_oscillator = false;
            d2bleWell = false;
            IMPT_PE=false;
            if app.ImportedPotentialCheckBox.Value
                IMPT_PE=true;
            end
            if app.FiniteSquareWellCheckBox.Value
                FiniteSquareWell = true;
            end
            if app.FiniteSquareBarrierCheckBox.Value
                SquareBarrier = true;
            end
            if app.CoulombPotentialCheckBox.Value
                CoulmbicInteractionPotential = true;
            end
            if app.StepPotentialCheckBox.Value
                StepFuncPotential = true;
            end
            if app.HarmonicOscillatorCheckBox.Value
                HarmonicOscillator = true;
            end
            if app.RandomPotentialCheckBox.Value
                RandomNoise = true;
            end
            if app.MorseOscillatorCheckBox.Value
                morse_oscillator = true;
            end
            if app.DoubleWellCheckBox.Value
                d2bleWell=true;
            end
            %% Kinetic Energy Matrix

            kinetic = zeros(DataPointsNumber);
            
            for i = 1:DataPointsNumber
                kinetic(i,i) = 2;
                if i > 1
                    kinetic(i, i-1) = -1;
                    kinetic(i-1, i) = -1;
                end
            end
            
            if PerBoundCond
                kinetic(DataPointsNumber,1) = -1;
                kinetic(1,DataPointsNumber) = -1;
            end
            
            kinetic_multiplier = hbar^2/(2*mass*step_spacing^2);
            kinetic = kinetic*kinetic_multiplier;
            
            %% Calculate Potential Energy Matrix
  
            
            potential = zeros(DataPointsNumber);
            
            
            %             app.CloumbicasymptotyicValue
            if CoulmbicInteractionPotential
                zz=app.FirstParticleCharge*app.SecondParticleCharge;
                for i = 1:DataPointsNumber
                    potential(i,i) = app.CloumbicasymptotyicValue+ (zz)/(i*step_spacing) + potential(i,i);
                end
            end
            
            if HarmonicOscillator
                Equ=round(app.EquPositionEditField.Value.*DataPointsNumber);
                FrcCnstnt=app.ForceConstantEditField.Value;
                Const=FrcCnstnt./2;
                for i = 1:DataPointsNumber
                    potential(i,i) = Const*(i - Equ)^2 + potential(i,i);
                end
            end
            
            
            % To consider MaximumAllowedCloumbic_Harmonic Value of potential Energy
            if CoulmbicInteractionPotential || HarmonicOscillator
                
                for i = 1:DataPointsNumber
                    if potential(i,i) > app.MaximumAllowedCloumbic_Harmonic
                        potential(i,i) = app.MaximumAllowedCloumbic_Harmonic;
                    end
                end
            end
            
            if IMPT_PE
                if app.IMP_TRUE
                    IMP_PTENGY1=app.IPF(:,2);
                    if DataPointsNumber~=1000
                        XIMP=app.IPF(:,1);
                        x=linspace(0,step_spacing*DataPointsNumber,DataPointsNumber);x=x';
                        IMP_PTENGY=spline(XIMP,IMP_PTENGY1,x);
                    else
                        IMP_PTENGY=IMP_PTENGY1;
                    end
                    for i = 1:DataPointsNumber
                        potential(i,i)=potential(i,i)+IMP_PTENGY(i,1);
                    end
                else
                    msgbox('No potential is imported or maybe the file format is not supported')
                end
            end
            
            
            
            if morse_oscillator
                de = app.DissociationEnergyEditField.Value;%0.1816; %H2 dissociation in hartrees
                k=app.ForceConstantEditField_2.Value; % k = 0.363; %force constant in hartrees
                alpha = sqrt(k/(2*de));
                Equ2=round(app.EquPositionEditField_2.Value.*DataPointsNumber);
                for i = 1:DataPointsNumber
                    x = 0.01*(i - Equ2);
                    potential(i,i) = de*((1-exp(-alpha*x))^2) + potential(i,i);
                    if potential(i,i) >= max(app.nNnTimes_morse*de,app.MaxLevelMorse)
                        potential(i,i) = max(app.nNnTimes_morse*de,app.MaxLevelMorse);
                    end
                end
            end
            
            if StepFuncPotential
                WP=round(app.WallPositionEditField.Value.*DataPointsNumber);
                RP=app.RightPotentialEditField.Value;
                LP=app.LeftPotentialEditField.Value;
                for i = 1:WP
                    potential(i,i) = LP + potential(i,i);
                end
                for i = WP + 1:DataPointsNumber
                    potential(i,i) = RP + potential(i,i);
                end
            end
            
            if FiniteSquareWell
                
                RW =round(app.RightWallPositionEditField.Value.*DataPointsNumber);
                LW =round(app.LeftWallPositionEditField.Value.*DataPointsNumber);
                Depth=app.WellDepthEditField.Value;
                for i = 1:LW
                    potential(i,i) = Depth + potential(i,i);
                end
                for i = LW+1:RW
                    potential(i,i) = 0 + potential(i,i);
                end
                for i = RW+1:DataPointsNumber
                    potential(i,i) = Depth + potential(i,i);
                end
            end
            
            
            if SquareBarrier
                LW1 =round(app.LeftWallPositionEditField_2.Value.*DataPointsNumber);
                RW1 =round(app.RightWallPositionEditField_2.Value.*DataPointsNumber);
                BH  =app.BarrierHeightEditField.Value;
                for i = 1:LW1
                    potential(i,i) = 0.0 + potential(i,i);
                end
                for i = LW1+1:RW1
                    potential(i,i) = BH + potential(i,i);
                end
                for i = RW1+1:DataPointsNumber
                    potential(i,i) = 0.0 + potential(i,i);
                end
            end
            
            if RandomNoise
                w=app.NoiseWeightEditField.Value;
                %                 potential(1,1) = 0;
                for i = 1:DataPointsNumber
                    potential(i,i) = w*(2*rand-1) + potential(i,i);
                    if potential(i,i) < 0
                        potential(i,i) = 0;
                    end
                end
                for i = 2:DataPointsNumber-1
                    potential(i,i) = (potential(i-1,i-1)+2*potential(i,i)+potential(i+1,i+1))/4;
                end
            end
            if d2bleWell
                x1r=app.x1EditField.Value; x1=round(x1r.*DataPointsNumber);
                x2r=app.x2EditField.Value; x2=round(x2r.*DataPointsNumber);
                x3r=app.x3EditField.Value; x3=round(x3r.*DataPointsNumber);
                x4r=app.x4EditField.Value; x4=round(x4r.*DataPointsNumber);
                a=app.aEditField.Value;
                c=app.cEditField.Value;
                PEM=zeros(size(potential,1),1);
                for i = 1:DataPointsNumber
                    PE=c.*1.0e-10.*(a.*i-x1).*(i-x2).*(i-x3).*(i-x4);
                    PEM(i)=PE;
                    potential(i,i) = PE  + potential(i,i);
                    
                end
                if app.setMinmum2Zero2Well
                    minval=min(PEM);
                    for i = 1:DataPointsNumber
                        potential(i,i)=  potential(i,i)-minval  ;
                    end
                end
            end
            if ~app.AllowNegatives
                for i=1:DataPointsNumber
                    if potential(i,i) < 0   %%% Should always be after the last potential
                        potential(i,i) = 0;
                    end
                end
            end
            %% Compute Eigenvectors of Hamiltonian

            app.Finalpotential=potential;
            hamiltonian = kinetic + potential;
            [eigenvectors, eigenvalues] = eig(hamiltonian);
            app.EigVal=eigenvalues;
            app.EigVec=eigenvectors;
            %% Plot Eigenvectors
             
            x_positions = linspace(0,step_spacing*DataPointsNumber,DataPointsNumber);
            
            %ground_state_wavefunction = eigenvectors(:,1) + eigenvalues(1,1)*ones(DataPointsNumber,1);
            %first_excited_state = eigenvectors(:,2) + eigenvalues(2,2)*ones(DataPointsNumber,1);
            %second_excited_state = eigenvectors(:,3) + eigenvalues(3,3)*ones(DataPointsNumber,1);
            %fifth_excited_state = eigenvectors(:,4) + eigenvalues(4,4)*ones(DataPointsNumber,1);
            if app.AddEigv2EigfPlot
                switch string(app.MultiSwitch.Value)
                    case "Off"
                        ChosenState=app.EigenStateQuantumNumberSpinner.Value;
                        if app.NormalizationSwitch
                            Normaa=trapz(x_positions,(eigenvectors(:,ChosenState)).^2);
                        else
                            Normaa=app.NormalizingConstant;
                        end
                        plot_wavefunction = eigenvectors(:,ChosenState)./sqrt(Normaa) + eigenvalues(ChosenState,ChosenState)*ones(DataPointsNumber,1);
                        area(app.UIAxes,x_positions, diag(potential),"DisplayName","Potential Energy")
                        hold(app.UIAxes,"on");
                        string2show=sprintf('# %i; Eigenvalue: %6.2f Hartrees',ChosenState,eigenvalues(ChosenState,ChosenState));
                        plot(app.UIAxes,x_positions, plot_wavefunction, 'g', 'LineWidth', 2,"DisplayName",string2show)
                        if app.Wavefunction2Switch.Value=="On"
                            plot_wavefunction2=((eigenvectors(:,ChosenState)).^2)./Normaa + eigenvalues(ChosenState,ChosenState)*ones(DataPointsNumber,1);
                            string2shows=sprintf('|Wavefunction# %i|%c',ChosenState,178);
                            plot(app.UIAxes,x_positions, plot_wavefunction2, 'k', 'LineWidth', 2,"DisplayName",string2shows)
                        end
                        xlabel(app.UIAxes,'Position')
                        ylabel(app.UIAxes,'Energy and Wavefunction')
                        %title('Wavefunction vs. Position for the First Three Energy Levels')
                        hold(app.UIAxes,"off");
                        %eigenvalues(1:20,1:20)
                        legend(app.UIAxes);
                    case "On"
                        state1=app.StateSpinner.Value;
                        state2=app.StateSpinner_2.Value;
                        state3=app.StateSpinner_3.Value;
                        state4=app.StateSpinner_4.Value;
                        states=[state1,state2,state3,state4];
                        if any(states)
                            area(app.UIAxes,x_positions, diag(potential),"DisplayName","Potential Energy")
                            xlabel(app.UIAxes,'Position')
                            ylabel(app.UIAxes,'Energy and Wavefunction')
                            
                        end
                        if state1>0
                            if app.NormalizationSwitch
                                Norm1=trapz(x_positions,(eigenvectors(:,state1)).^2);
                            else
                                Norm1=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state1)./sqrt(Norm1) + eigenvalues(state1,state1)*ones(DataPointsNumber,1);
                            hold(app.UIAxes,"on");
                            string2show1=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state1,eigenvalues(state1,state1));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'k', 'LineWidth', 2,"DisplayName",string2show1)
                        end
                        if state2>0
                            if app.NormalizationSwitch
                                Norm2=trapz(x_positions,(eigenvectors(:,state2)).^2);
                            else
                                Norm2=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state2)./sqrt(Norm2) + eigenvalues(state2,state2)*ones(DataPointsNumber,1);
                            hold(app.UIAxes,"on");
                            string2show2=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state2,eigenvalues(state2,state2));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'r', 'LineWidth', 2,"DisplayName",string2show2)
                        end
                        if state3>0
                            if app.NormalizationSwitch
                                Norm3=trapz(x_positions,(eigenvectors(:,state3)).^2);
                            else
                                Norm3=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state3)./sqrt(Norm3) + eigenvalues(state3,state3)*ones(DataPointsNumber,1);
                            hold(app.UIAxes,"on");
                            string2show3=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state3,eigenvalues(state3,state3));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'b', 'LineWidth', 2,"DisplayName",string2show3)
                        end
                        if state4>0
                            if app.NormalizationSwitch
                                Norm4=trapz(x_positions,(eigenvectors(:,state4)).^2);
                            else
                                Norm4=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state4)./sqrt(Norm4) + eigenvalues(state4,state4)*ones(DataPointsNumber,1);
                            hold(app.UIAxes,"on");
                            string2show4=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state4,eigenvalues(state4,state4));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'g', 'LineWidth', 2,"DisplayName",string2show4)
                        end
                        legend(app.UIAxes);
                        hold(app.UIAxes,"off");
                        
                end
                
            else
                switch string(app.MultiSwitch.Value)
                    case "Off"
                        ChosenState=app.EigenStateQuantumNumberSpinner.Value;
                        if app.NormalizationSwitch
                            Normaa=trapz(x_positions,(eigenvectors(:,ChosenState)).^2);
                        else
                            Normaa=app.NormalizingConstant;
                        end
                        plot_wavefunction = eigenvectors(:,ChosenState)./sqrt(Normaa) ;
                        area(app.UIAxes,x_positions, diag(potential),"DisplayName","Potential Energy")
                        hold(app.UIAxes,"on");
                        string2show=sprintf('# %i; Eigenvalue: %6.2f',ChosenState,eigenvalues(ChosenState,ChosenState));
                        plot(app.UIAxes,x_positions, plot_wavefunction, 'g', 'LineWidth', 2,"DisplayName",string2show)
                        if app.Wavefunction2Switch.Value=="On"
                            plot_wavefunction2=(eigenvectors(:,ChosenState).^2)./Normaa;
                            string2shows=sprintf('|Wavefunction# %i|%c',ChosenState,178);
                            plot(app.UIAxes,x_positions, plot_wavefunction2, 'k', 'LineWidth', 2,"DisplayName",string2shows)
                        end
                        xlabel(app.UIAxes,'Position')
                        ylabel(app.UIAxes,'Energy and Wavefunction')
                        %title('Wavefunction vs. Position for the First Three Energy Levels')
                        hold(app.UIAxes,"off");
                        %eigenvalues(1:20,1:20)
                        legend(app.UIAxes);
                    case "On"
                        state1=app.StateSpinner.Value;
                        state2=app.StateSpinner_2.Value;
                        state3=app.StateSpinner_3.Value;
                        state4=app.StateSpinner_4.Value;
                        states=[state1,state2,state3,state4];
                        if any(states)
                            area(app.UIAxes,x_positions, diag(potential),"DisplayName","Potential Energy")
                            xlabel(app.UIAxes,'Position')
                            ylabel(app.UIAxes,'Energy and Wavefunction')
                            
                        end
                        if state1>0
                            if app.NormalizationSwitch
                                Norm1=trapz(x_positions,(eigenvectors(:,state1)).^2);
                            else
                                Norm1=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state1)./sqrt(Norm1) ;
                            hold(app.UIAxes,"on");
                            string2show1=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state1,eigenvalues(state1,state1));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'k', 'LineWidth', 2,"DisplayName",string2show1)
                        end
                        if state2>0
                            if app.NormalizationSwitch
                                Norm2=trapz(x_positions,(eigenvectors(:,state2)).^2);
                            else
                                Norm2=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state2)./sqrt(Norm2);
                            hold(app.UIAxes,"on");
                            string2show2=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state2,eigenvalues(state2,state2));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'r', 'LineWidth', 2,"DisplayName",string2show2)
                        end
                        if state3>0
                            if app.NormalizationSwitch
                                Norm3=trapz(x_positions,(eigenvectors(:,state3)).^2);
                            else
                                Norm3=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state3)./sqrt(Norm3);
                            hold(app.UIAxes,"on");
                            string2show3=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state3,eigenvalues(state3,state3));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'b', 'LineWidth', 2,"DisplayName",string2show3)
                        end
                        if state4>0
                            if app.NormalizationSwitch
                                Norm4=trapz(x_positions,(eigenvectors(:,state4)).^2);
                            else
                                Norm4=app.NormalizingConstant;
                            end
                            plot_wavefunction1 = eigenvectors(:,state4)./sqrt(Norm4) ;
                            hold(app.UIAxes,"on");
                            string2show4=sprintf('# %i; Eigenvalue: %6.2f Hartrees',state4,eigenvalues(state4,state4));
                            plot(app.UIAxes,x_positions, plot_wavefunction1, 'g', 'LineWidth', 2,"DisplayName",string2show4)
                        end
                        legend(app.UIAxes);
                        hold(app.UIAxes,"off");
                        
                end
                
            end
        end

        % Image clicked function: Image
        function ImageClicked(app, event)
            msgbox('Nima Soltani and Professor John S. Hutchinson; Rice University; 2021');
            if strcmp(app.Image4.Visible,'off')
                app.Image4.Visible='on';
            else
                app.Image4.Visible='off';
            end
        end

        % Value changed function: MultiSwitch
        function MultiSwitchValueChanged(app, event)
            if strcmp(string(app.MultiSwitch.Value),"Off")
                app.StateSpinner.BackgroundColor='k';
                app.StateSpinner_2.BackgroundColor='r';
                app.StateSpinner_3.BackgroundColor='b';
                app.StateSpinner_4.BackgroundColor='g';
                app.Wavefunction2Switch.Visible='on';
                app.Wavefunction2Switch.Enable='on';
                app.MultiStateVisualizationmodeLabel.BackgroundColor=[0.94,0.94,0.94];
                app.MultiStateVisualizationmodeLabel.Visible='off';
                app.EigenStateQuantumNumberSpinner.BackgroundColor=[1.00,1.0,1.0];
            elseif strcmp(string(app.MultiSwitch.Value),"On")
                app.StateSpinner.BackgroundColor=[1.00,1.0,1.0];
                app.StateSpinner_2.BackgroundColor=[1.00,1.0,1.0];
                app.StateSpinner_3.BackgroundColor=[1.00,1.0,1.0];
                app.StateSpinner_4.BackgroundColor=[1.00,1.0,1.0];
                app.EigenStateQuantumNumberSpinner.BackgroundColor='k';
                app.Wavefunction2Switch.Visible='off';
                app.MultiStateVisualizationmodeLabel.Visible='on';
                app.Wavefunction2Switch.Enable='off';
                app.MultiStateVisualizationmodeLabel.BackgroundColor='y';
            end
        end

        % Button pushed function: ResettoDefaultButton
        function ResettoDefaultButtonPushed(app, event)
            app.LeftWallPositionEditField.Limits=[0,1];
            app.RightWallPositionEditField.Limits=[0,1];
            app.LeftWallPositionEditField.Value=0.4;
            app.RightWallPositionEditField.Value=0.6;
            app.WellDepthEditField.Value=0.9;
        end

        % Button pushed function: ResettoDefaultButton_2
        function ResettoDefaultButton_2Pushed(app, event)
            app.LeftWallPositionEditField_2.Limits=[0,1];
            app.RightWallPositionEditField_2.Limits=[0,1];
            app.LeftWallPositionEditField_2.Value=0.4;
            app.RightWallPositionEditField_2.Value=0.6;
            app.BarrierHeightEditField.Value=1.0;
        end

        % Value changed function: FiniteSquareBarrierCheckBox
        function FiniteSquareBarrierCheckBoxValueChanged(app, event)
            if app.FiniteSquareBarrierCheckBox.Value
                app.LeftWallPositionEditField_2.Enable='on';
                app.RightWallPositionEditField_2.Enable='on';
                app.BarrierHeightEditField.Enable='on';
                app.ResettoDefaultButton_2.Enable='on';
            else
                app.LeftWallPositionEditField_2.Enable='off';
                app.RightWallPositionEditField_2.Enable='off';
                app.BarrierHeightEditField.Enable='off';
                app.ResettoDefaultButton_2.Enable='off';
            end
        end

        % Button pushed function: ResettoDefaultButton_3
        function ResettoDefaultButton_3Pushed(app, event)
            app.WallPositionEditField.Value=0.5;
            app.RightPotentialEditField.Value=1;
            app.LeftPotentialEditField.Value=0.1;
        end

        % Value changed function: StepPotentialCheckBox
        function StepPotentialCheckBoxValueChanged(app, event)
            if app.StepPotentialCheckBox.Value
                app.WallPositionEditField.Enable='on';
                app.RightPotentialEditField.Enable='on';
                app.LeftPotentialEditField.Enable='on';
                app.ResettoDefaultButton_3.Enable='on';
            else
                app.WallPositionEditField.Enable='off';
                app.RightPotentialEditField.Enable='off';
                app.LeftPotentialEditField.Enable='off';
                app.ResettoDefaultButton_3.Enable='off';
            end
        end

        % Value changed function: HarmonicOscillatorCheckBox
        function HarmonicOscillatorCheckBoxValueChanged(app, event)
            if app.HarmonicOscillatorCheckBox.Value
                app.EquPositionEditField.Enable='on';
                app.ForceConstantEditField.Enable='on';
                app.ResettoDefaultButton_4.Enable='on';
            else
                app.EquPositionEditField.Enable='off';
                app.ForceConstantEditField.Enable='off';
                app.ResettoDefaultButton_4.Enable='off';
            end
            
        end

        % Value changed function: MorseOscillatorCheckBox
        function MorseOscillatorCheckBoxValueChanged(app, event)
            if app.MorseOscillatorCheckBox.Value
                app.ResettoDefaultButton_5.Enable='on';
                app.ForceConstantEditField_2.Enable='on';
                app.EquPositionEditField_2.Enable='on';
                app.DissociationEnergyEditField.Enable='on';
            else
                app.ResettoDefaultButton_5.Enable='off';
                app.ForceConstantEditField_2.Enable='off';
                app.EquPositionEditField_2.Enable='off';
                app.DissociationEnergyEditField.Enable='off';
            end
        end

        % Button pushed function: ResettoDefaultButton_5
        function ResettoDefaultButton_5Pushed(app, event)
            app.ForceConstantEditField_2.Value=0.363;
            app.EquPositionEditField_2.Value=0.5;
            app.DissociationEnergyEditField.Value=0.1816;
        end

        % Value changed function: DoubleWellCheckBox
        function DoubleWellCheckBoxValueChanged(app, event)
            if app.DoubleWellCheckBox.Value
                app.ResettoDefaultButton_6.Enable='on';
                app.aEditField.Enable='on';
                app.x1EditField.Enable='on';
                app.x2EditField.Enable='on';
                app.x3EditField.Enable='on';
                app.x4EditField.Enable='on';
                app.cEditField.Enable='on';
                app.PEcaxx1xx2xx3xx4Label.Enable='on';
            else
                app.ResettoDefaultButton_6.Enable='off';
                app.aEditField.Enable='off';
                app.x1EditField.Enable='off';
                app.x2EditField.Enable='off';
                app.x3EditField.Enable='off';
                app.x4EditField.Enable='off';
                app.PEcaxx1xx2xx3xx4Label.Enable='off';
                app.cEditField.Enable='off';
            end
        end

        % Button pushed function: ResettoDefaultButton_6
        function ResettoDefaultButton_6Pushed(app, event)
            app.aEditField.Value=1;
            app.x1EditField.Value=0.02;
            app.x2EditField.Value=0.35;
            app.x3EditField.Value=0.65;
            app.x4EditField.Value=0.98;
            app.cEditField.Value=1;
        end

        % Menu selected function: ShowmemoreMenu
        function ShowmemoreMenuSelected(app, event)
            app.Image3.Visible='on';
            app.Image2.Visible='on';
        end

        % Menu selected function: FlythemoffMenu
        function FlythemoffMenuSelected(app, event)
            app.Image4.Visible='off';
            app.Image3.Visible='off';
            app.Image2.Visible='off';
        end

        % Menu selected function: DefaultMenu
        function DefaultMenuSelected(app, event)
            app.Image4.Visible='off';
            app.Image3.Visible='off';
            app.Image2.Visible='off';
            app.UIFigure.Color=[0.94,0.94,0.94];
            app.UIAxes.BackgroundColor=[0.94,0.94,0.94];
            %             app.RunButton.BackgroundColor=[1,1,1];
        end

        % Menu selected function: BlueMenu
        function BlueMenuSelected(app, event)
            app.UIFigure.Color=[0.30,0.75,0.93];
            app.UIAxes.BackgroundColor=[0.30,0.75,0.92];
            %             app.RunButton.BackgroundColor=[0.30,0.75,0.93];
        end

        % Menu selected function: WhiteMenu
        function WhiteMenuSelected(app, event)
            app.UIFigure.Color=[1,1,1];
            app.RunButton.BackgroundColor=[1,1,1];
            app.UIAxes.BackgroundColor=[1,1,1];
        end

        % Menu selected function: DarkMenu
        function DarkMenuSelected(app, event)
            app.UIFigure.Color=[0.42,0.42,0.42];
            app.UIAxes.BackgroundColor=[0.45,0.45,0.45];
            %             app.RunButton.BackgroundColor=[0.7,0.7,0.7];
        end

        % Callback function
        function JPGMenuSelected(app, event)
            
        end

        % Button pushed function: ResettoDefaultButton_7
        function ResettoDefaultButton_7Pushed(app, event)
            app.NoiseWeightEditField.Value=0.1;
        end

        % Value changed function: RandomPotentialCheckBox
        function RandomPotentialCheckBoxValueChanged(app, event)
            if app.RandomPotentialCheckBox.Value
                app.NoiseWeightEditField.Enable='on';
                app.ResettoDefaultButton_7.Enable='on';
            else
                app.ResettoDefaultButton_7.Enable='off';
                app.NoiseWeightEditField.Enable='off';
            end
            
        end

        % Value changed function: ImportedPotentialCheckBox
        function ImportedPotentialCheckBoxValueChanged(app, event)
            if app.ImportedPotentialCheckBox.Value
                app.ImportPotentialButton.Enable='on';
                app.IMPRTPETextArea.Enable='on';
                app.IMPRT_PELamp.Enable='on';
                if app.IMP_TRUE
                    app.IMPRT_PELamp.Color='g';
                    app.UIAxes2.Visible='on';
                else
                    app.IMPRT_PELamp.Color='r';
                    app.UIAxes2.Visible='off';
                end
            else
                app.ImportPotentialButton.Enable='off';
                app.IMPRTPETextArea.Enable='off';
                app.IMPRT_PELamp.Enable='off';
                app.IMPRT_PELamp.Color='b';
                app.UIAxes2.Visible='off';
            end
            
        end

        % Value changed function: ImportPotentialButton
        function ImportPotentialButtonValueChanged(app, event)
            % value = app.ImportPotentialButton.Value;
            prevpath=app.PRVpath;
            [app.IMPPE_File,app.IMPPE_Path,index] = uigetfile({'*.csv;*.xls;*.xlsb;*.xlsm;*.xlsx;*.xltm;*.xltx;*.ods','SPreadShit';'*.txt;*.csv;*.dat','Text & CSV Files';'*.*','All Files (*.*)'},'Import the Potential',prevpath);
            
            if index
                app.PRVpath=app.IMPPE_Path;
                filename1=[app.IMPPE_Path,app.IMPPE_File];
                app.IMPRTPETextArea.Value=filename1;
                readpotential=readmatrix(filename1);
                %app.IPF=csvread(filename1,0,0,[0 0 999 1]);
                SZ_C=size(readpotential,2);
                if SZ_C>0
                    xs=transpose(linspace(0,50,1000));
                    if SZ_C>2
                        readpotential=readpotential(:,1:2);
                        if size(readpotential,1)==1000
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                            msgbox(['Your data contains ',num2str(SZ_C),' columns; Only first two columns will be used.']);
                        elseif size(readpotential,1)>1000
                            msgbox(['Your data contains ',num2str(SZ_C),' columns; Only first two columns will be used. ','More than 1000 datapoints is found; first 1000 rows will be mapped tp 0-50']);
                            readpotential=readpotential(1:1000,1:2);
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        elseif size(readpotential,1)<5
                            msgbox('Your data has less than five data points, make sure your data is saved in two columns: Xs;Ys.  X ranging from 0 to 50');
                        else % which means it is between 5 datapoints and 1000 datapoints
                            msgbox(['Your data contains ',num2str(SZ_C),' columns; Only first two columns will be used. Less than 1000 datapoints is found; missing values will be set to zero']);
                            sz=size(readpotential,1);
                            cc=zeros([1000,2]);
                            cc(:,1)=xs;
                            cc(1:sz,2)=readpotential(:,2);
                            readpotential=cc;
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        end
                    elseif SZ_C==1
                        if size(readpotential,1)==1000
                            readpotential=[xs,readpotential];
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                            msgbox('It seems that your data has only one column, It will be used as Potential energy and it will be mapped to 0=<X=<50')
                        elseif size(readpotential,1)>1000
                            msgbox('More than 1000 datapoints is found; first 1000 rows will be mapped tp 0-50; It seems that your data has only one column, It will be used as Potential energy and it will be mapped to 0=<X=<50');
                            readpotential=readpotential(1:1000,1);
                            readpotential=[xs,readpotential];
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        elseif size(readpotential,1)<5
                            msgbox('Your data has less than five data points, make sure your data is saved in two columns: Xs;Ys.  X ranging from 0 to 50');
                        else % which means it is between 5 datapoints and 1000 datapoints
                            msgbox('Less than 1000 datapoints is found; missing values will be set to zero. It seems that your data has only one column, It will be used as Potential energy and it will be mapped to 0=<X=<50');
                            sz=size(readpotential,1);
                            cc=zeros([1000,2]);
                            cc(:,1)=xs;
                            cc(1:sz,2)=readpotential;
                            readpotential=cc;
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        end
                        
                    elseif SZ_C==2
                        if size(readpotential,1)==1000
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        elseif size(readpotential,1)>1000
                            msgbox('More than 1000 datapoints is found; first 1000 rows will be mapped tp 0-50');
                            readpotential=readpotential(1:1000,1:2);
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        elseif size(readpotential,1)<5
                            msgbox('Your data has less than five data points, make sure your data is saved in two columns: Xs;Ys.  X ranging from 0 to 50');
                        else % which means it is between 5 datapoints and 1000 datapoints
                            msgbox('Less than 1000 datapoints is found; missing values will be set to zero');
                            sz=size(readpotential,1);
                            cc=zeros([1000,2]);
                            cc(:,1)=xs;
                            cc(1:sz,2)=readpotential(:,2);
                            readpotential=cc;
                            app.IPF=readpotential;
                            app.IMPRT_PELamp.Color='g';
                            app.IMP_TRUE=1;
                            app.UIAxes2.Visible='on';
                            area(app.UIAxes2,xs,readpotential(:,2))
                        end
                    end
                    
                else
                    msgbox('No Potential Is Imported;  make sure your data is saved in two columns with 1000 data points (Xs,Ys); X ranging from 0 to 50');
                end
            end
        end

        % Menu selected function: ExportPotentialMenu_3
        function ExportPotentialMenu_3Selected(app, event)
            prevpath=app.PRVpath;
            [app.ExportFileName2,app.ExportPathName2,index]=uiputfile({'*.csv','.CSV';'*.txt;*.csv','TXT & CSV Files';'*.*','All Files (*.*)'},'Export Potential',prevpath);
            if index
                step_spacing = app.StepSize;
                app.PRVpath=app.ExportPathName2;
                DataPointsNumber = app.Npoints;
                if DataPointsNumber==1000
                    PE2Write=diag(app.Finalpotential);
                    x=linspace(0,step_spacing*DataPointsNumber,DataPointsNumber);
                    X2Write=x';
                else
                    PE2extrapolate=diag(app.Finalpotential);
                    x=linspace(0,step_spacing*DataPointsNumber,DataPointsNumber);
                    x1=linspace(0,0.05*1000,1000);
                    PE2Write=spline(x,PE2extrapolate,x1);PE2Write=PE2Write';
                    X2Write=x1';
                end
                EXPORT_PE=[X2Write,PE2Write];
                filename=[app.ExportPathName2,app.ExportFileName2];
                writematrix(EXPORT_PE,filename);
            end
            
        end

        % Value changed function: Button
        function ButtonValueChanged(app, event)
            %             value = app.Button.Value;
            msgbox('2 column by 1000 row (Xs,Ys) Data: Import data in the format of spreadsheet (.xls, .xlsb, .xlsm, .xlsx, .xltm, .xltx, or .ods) or  delimited text files (.txt, .dat, or .csv ), X ranging from 0 to 50-- with a separation of 50/999. Note if the data contains 1 column, it is interpreted as Y data and it will be mapped to x 0-50; in fact every sets of acceptable data (2 coulmn by 1000 rows) is also mapped to x [0-50] regardless of their x values!')
            
        end

        % Menu selected function: AboutMenu
        function AboutMenuSelected(app, event)
            msgbox('GNU Liscence. This app is developed only for educational purposes and it is not meant to be sold or purchased. For suggestions please contact Nima Soltani via nima.slt@gmail.com or Professor John S. Hutchinson at Rice university');
        end

        % Menu selected function: DocumentationMenu
        function DocumentationMenuSelected(app, event)
            msgbox('For more information please refer to https://sites.google.com/view/schrodingerwave ; Units are in Hartrees. Changing resolution will affect speed. Program is designed to enable sharing Potential--comsistent between differnet resolutions. ');
        end

        % Menu selected function: SettoHighMenu
        function SettoHighMenuSelected(app, event)
            app.Npoints=5000;
            app.StepSize=0.01;
            app.EigenStateQuantumNumberSpinner.Limits=[1,app.Npoints];
            app.StateSpinner.Limits=[0,app.Npoints];
            app.StateSpinner_2.Limits=[0,app.Npoints];
            app.StateSpinner_3.Limits=[0,app.Npoints];
            app.StateSpinner_4.Limits=[0,app.Npoints];
            msgbox('Resolution Set to High. Unites are in Hartree')
        end

        % Menu selected function: SettoMediumMenu
        function SettoMediumMenuSelected(app, event)
            app.Npoints=1000;
            app.StepSize=5/100;
            app.EigenStateQuantumNumberSpinner.Limits=[1,app.Npoints];
            app.StateSpinner.Limits=[0,app.Npoints];
            app.StateSpinner_2.Limits=[0,app.Npoints];
            app.StateSpinner_3.Limits=[0,app.Npoints];
            app.StateSpinner_4.Limits=[0,app.Npoints];
            msgbox('Resolution Set to Medium. Unites are in Hartree')
        end

        % Menu selected function: SettoLowMenu
        function SettoLowMenuSelected(app, event)
            app.Npoints=200;
            app.StepSize=5/20;
            app.EigenStateQuantumNumberSpinner.Limits=[1,app.Npoints];
            app.StateSpinner.Limits=[0,app.Npoints];
            app.StateSpinner_2.Limits=[0,app.Npoints];
            app.StateSpinner_3.Limits=[0,app.Npoints];
            app.StateSpinner_4.Limits=[0,app.Npoints];
            msgbox('Resolution Set to Low. Unites are in Hartree')
        end

        % Value changed function: LeftWallPositionEditField
        function LeftWallPositionEditFieldValueChanged(app, event)
            LWell = app.LeftWallPositionEditField.Value;
            app.RightWallPositionEditField.Limits=[LWell,1];
        end

        % Value changed function: RightWallPositionEditField
        function RightWallPositionEditFieldValueChanged(app, event)
            RWell = app.RightWallPositionEditField.Value;
            app.LeftWallPositionEditField.Limits=[0,RWell];
        end

        % Value changed function: LeftWallPositionEditField_2
        function LeftWallPositionEditField_2ValueChanged(app, event)
            LBarrier = app.LeftWallPositionEditField_2.Value;
            app.RightWallPositionEditField_2.Limits=[LBarrier,1];
        end

        % Value changed function: RightWallPositionEditField_2
        function RightWallPositionEditField_2ValueChanged(app, event)
            RBarrier = app.RightWallPositionEditField_2.Value;
            app.LeftWallPositionEditField_2.Limits=[0,RBarrier];
        end

        % Value changed function: RightPotentialEditField
        function RightPotentialEditFieldValueChanged(app, event)
            %value = app.RightPotentialEditField.Value;
            % redundant
        end

        % Menu selected function: MorseOscillatorMenu
        function MorseOscillatorMenuSelected(app, event)
            % redundant
        end

        % Menu selected function: SettingMenu_2
        function SettingMenu_2Selected(app, event)
            %       nNnTimes_morse =2.5;        MaxLevelMorse=0.5;
            definput1 = {num2str(app.nNnTimes_morse),num2str(app.MaxLevelMorse)};
            prompt = {'N multiple of Dissociation Energy','Preset maximum value of Energy (Hartree)'};
            dlgtitle = 'The Maximum Allowed Energy is the Maximum of these values:';
            dims = [1 70];
            answer = inputdlg(prompt,dlgtitle,dims,definput1);
            if max(size(answer))>=1
                app.nNnTimes_morse=str2double(answer{1,1});
                app.MaxLevelMorse=str2double(answer{2,1});
            end
        end

        % Menu selected function: HelpMenu_2
        function HelpMenu_2Selected(app, event)
            msgbox('Since energy of Morse Oscillator blows up for low values of x, we put a maximum allowed value for energy. The maximum value of allowed energy is either a preset level "Preset maximum value of Energy (Hartree)" or multiple of De: "N multiple of Dissociation Energy", whichever is  greater');
        end

        % Menu selected function: ShiftWavefunctionbyEnergyMenu
        function ShiftWavefunctionbyEnergyMenuSelected(app, event)
            %%% Everuthing transfered to  EnableorDisableMenuSelected(app, event)
            
            %             Listf={'Add Eigenvalue','Do not add Eigenvalue'};
            %             [Index,~]=listdlg('PromptString',{'Select a mode for plotting.',...
            %              'Vertically Shifting Eigenfunction by Adding Corresponding Eigenvalue or Not','','',''}, 'SelectionMode','single','ListString',Listf);
            %          if Index==1
            %              app.AddEigv2EigfPlot=1;
            %          elseif Index==2
            %              app.AddEigv2EigfPlot=0;
            %          end
            %
        end

        % Menu selected function: SelectandPlotEigenvaluesMenu
        function SelectandPlotEigenvaluesMenuSelected(app, event)
            ss=1;
            nn=25;
            definput1 = {num2str(ss),num2str(nn)};
            prompt = {'Lower limit State# (minimum =1)',['Upper limit State# (maximum =',num2str(app.Npoints),')']};
            dlgtitle = 'Range of Eigenstates to plot energy';
            dims = [1 70];
            answer = inputdlg(prompt,dlgtitle,dims,definput1);
            if max(size(answer))>=1
                ss=ceil(str2double(answer{1,1}));
                nn=floor(str2double(answer{2,1}));
                if (ss>0 && ss<=app.Npoints) && (nn>=ss && nn<=app.Npoints)
                    ffig = uifigure;
                    ff=uiaxes(ffig);
                    eigenvalues=app.EigVal;
                    ddd = diag(eigenvalues);
                    First25=ddd(ss:nn);
                    x=ss:1:nn;
                    %                 x=ones(size(First25));
                    %                 for ii=ss:nn
                    %                     x(ii)=ii;
                    %                 end
                    %
                    if app.PlotType<=2
                        b=bar(ff,x,First25);
                        if app.PlotType==1
                            xtips1 = b.XEndPoints;
                            ytips1 = b.YEndPoints;
                            labels1 = string(b.YData);
                            text(ff,xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
                        end
                    elseif app.PlotType==3
                        plot(ff,x,First25,'b-');
                    elseif app.PlotType==4
                        plot(ff,x,First25,'bo');
                    elseif app.PlotType==5
                        plot(ff,x,First25,'bo-');
                    end
                    ylabel(ff,'Energy (Hartrees)')
                    xlabel(ff,'State #')
                else
                    msgbox('Please choose your range within the limits')
                end
            end
        end

        % Menu selected function: EnableorDisableMenu
        function EnableorDisableMenuSelected(app, event)
            Listf={'Add Eigenvalue','Do not add Eigenvalue'};
            [Index,~]=listdlg('PromptString',{'Select a mode for plotting.',...
                'Vertically Shifting Eigenfunction by Adding Corresponding Eigenvalue or Not','','',''}, 'SelectionMode','single','ListString',Listf);
            if Index==1
                app.AddEigv2EigfPlot=1;
            elseif Index==2
                app.AddEigv2EigfPlot=0;
            end
        end

        % Menu selected function: ExportEigenvectorscsvMenu
        function ExportEigenvectorscsvMenuSelected(app, event)
            eigenvectors=app.EigVec;
            prevpath=app.PRVpath;
            [file,path,indx] = uiputfile('*.csv','Save Eigenvectors',prevpath);
            if indx
            app.PRVpath=path;
            filename=[path,file];
            writematrix(eigenvectors,filename)
            end
        end

        % Menu selected function: ExportEigenvaluescsvMenu_2
        function ExportEigenvaluescsvMenu_2Selected(app, event)
            eigenvals=diag(app.EigVal);
            prevpath=app.PRVpath;
            [file,path,indx] = uiputfile('*.csv','Save Eigenvalues',prevpath);
            if indx
            app.PRVpath=path;
            filename=[path,file];
            writematrix(eigenvals,filename)
            end
        end

        % Value changed function: Wavefunction2Switch
        function Wavefunction2SwitchValueChanged(app, event)
            %             value = app.Wavefunction2Switch.Value;
            
        end

        % Menu selected function: ChangeNormalizationConstantMenu
        function ChangeNormalizationConstantMenuSelected(app, event)
            Listf={'Automatic Normalization','Do not Normalize'};
            [Index,~]=listdlg('PromptString',{'Select a mode:','',...
                'To  defined the value of a constant multiplier for all the wavefunctions, select "Do not Normalize" ','','',''}, 'SelectionMode','single','ListString',Listf);
            if Index==1
                app.NormalizationSwitch=1;
                %                  app.NormalizingConstant=1;
            elseif Index==2
                app.NormalizationSwitch=0;
                
            end
            if app.NormalizationSwitch==0
                definput1 = {num2str(sqrt(1./app.NormalizingConstant))};
                prompt = {'All Wavefunctions will be multiplied to this constants (Or its squared value in the case of squared wavefunction)'};
                dlgtitle = 'Define Constant';
                dims = [1 70];
                answer = inputdlg(prompt,dlgtitle,dims,definput1);
                if max(size(answer))>=1
                    aa=str2double(answer{1,1});
                    app.NormalizingConstant=1./(aa.^2);
                    msgbox('Auto Normalization Mode is Disabled')
                end
            else
                msgbox('Auto Normalization Mode is Enabled')
            end
        end

        % Menu selected function: SettingMenu_4
        function SettingMenu_4Selected(app, event)
            definput1 = {num2str(app.CloumbicasymptotyicValue),num2str(app.FirstParticleCharge),num2str(app.SecondParticleCharge)};
            prompt = {'The Asymptotic Value of potential energy when the charged particles are far apart (Hartrees)', 'Charge on the first particle (multiple of electron charge)','Charge on the second particle (multiple of electron charge)'};
            dlgtitle = 'Modify Cloumbic Interactions parameters';
            dims = [1 70];
            answer = inputdlg(prompt,dlgtitle,dims,definput1);
            if max(size(answer))>=1
                if str2double(answer{1,1})>app.MaximumAllowedCloumbic_Harmonic
                    msgbox('The Asymptotic Value of potential cannot be greater than the "Maximum Allowed Value of Potential Energy" for Columbic Interactions, You can change this value in the "Hidden setting of Harmonic oscillator"')
                else
                    app.CloumbicasymptotyicValue=str2double(answer{1,1});
                end
                app.FirstParticleCharge=str2double(answer{2,1});
                app.SecondParticleCharge=str2double(answer{3,1});
            end
        end

        % Menu selected function: HelpMenu_4
        function HelpMenu_4Selected(app, event)
            msgbox('By using the Columbic Hidden Setting menu, it is possible to change the charge of two interacting particles as well as the value of potential energy when the particles are far apart. In the case of repulsive columbic interactions, the maximum value of allowed energy (to prevent the energy blow up) can be adjusted through "{Harmonic Oscillator} Hidden Potental Setting". Note that Negative values for energy are not allowed.')
        end

        % Menu selected function: SettingMenu_3
        function SettingMenu_3Selected(app, event)
            definput1 = {num2str(app.MaximumAllowedCloumbic_Harmonic)};
            prompt = {'The maximum Allowed Value of Energy for Harmonic Oscillator and/or Columbic Interacations'};
            dlgtitle = 'Define Maximum Energy to prevent Blow up';
            dims = [1 70];
            answer = inputdlg(prompt,dlgtitle,dims,definput1);
            if max(size(answer))>=1
                app.MaximumAllowedCloumbic_Harmonic=str2double(answer{1,1});
            end
        end

        % Menu selected function: HelpMenu_3
        function HelpMenu_3Selected(app, event)
            msgbox('Harmonic Oscillator Hidden Potential Setting menu allows you to change the maximum allowed value of potential energy for either Cloumbic or Harmonic Oscillator. This setting will work if either of the H.O. or Cloumbic potentials is selected.  ')
        end

        % Menu selected function: InfoonResolutionMenu
        function InfoonResolutionMenuSelected(app, event)
            msgbox('High Res Mode= 5000 basis sets; Medium Res Mode=1000 basis sets; Low Res Mode=200 basis sets. The defulat setting is Medium. If you want to design a potential (from an extrnal app like Excel or MATLAB) Always do 1000 points, from 0 to 50-- with a separation of 50/999; the program can take care of converting it to match the resolution you have. Refer to the website for more information: https://sites.google.com/view/schrodingerwave ')
        end

        % Menu selected function: 
        % AllowNegativeValuesofPotentialEnergyMenu
        function AllowNegativeValuesofPotentialEnergyMenuSelected(app, event)
            definput1 = {num2str(app.AllowNegatives)};
            prompt = {'Allow Negative Values for Potential Energy? ("0" Zero=Does not allow;  "1" One =allow)'};
            dlgtitle = 'Negative Potential Energy?';
            dims = [1 70];
            answer = inputdlg(prompt,dlgtitle,dims,definput1);
            if max(size(answer))>=1
                app.AllowNegatives=str2double(answer{1,1});
                if app.AllowNegatives==0
                    msgbox('Negative values for potential energy will be replaced by zero; PE>0')
                else
                    msgbox('Potential Energy is allowed to assume negative values')
                end
            end
            
        end

        % Menu selected function: ChangePlotTypeMenu
        function ChangePlotTypeMenuSelected(app, event)
            Listf={'Bar (with labels)','Bar (without labels)','Line','symbols','Line and Symbols'};
            [Index,~]=listdlg('PromptString',{'Select a mode for plotting.',...
                'Tpye of Energy Spectrum Plot','','',''}, 'SelectionMode','single','ListString',Listf);
            if Index==1
                app.PlotType=1;
            elseif Index==2
                app.PlotType=2;
            elseif Index==3
                app.PlotType=3;
            elseif Index==4
                app.PlotType=4;
            elseif Index==5
                app.PlotType=5;
            else
                % Nothing :-)
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1066 676];
            app.UIFigure.Name = 'UI Figure';

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create ExportPotentialMenu_3
            app.ExportPotentialMenu_3 = uimenu(app.FileMenu);
            app.ExportPotentialMenu_3.MenuSelectedFcn = createCallbackFcn(app, @ExportPotentialMenu_3Selected, true);
            app.ExportPotentialMenu_3.Text = 'Export Potential';

            % Create ExportEigenvaluescsvMenu_2
            app.ExportEigenvaluescsvMenu_2 = uimenu(app.FileMenu);
            app.ExportEigenvaluescsvMenu_2.MenuSelectedFcn = createCallbackFcn(app, @ExportEigenvaluescsvMenu_2Selected, true);
            app.ExportEigenvaluescsvMenu_2.Text = 'Export Eigenvalues (.csv)';

            % Create ExportEigenvectorscsvMenu
            app.ExportEigenvectorscsvMenu = uimenu(app.FileMenu);
            app.ExportEigenvectorscsvMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportEigenvectorscsvMenuSelected, true);
            app.ExportEigenvectorscsvMenu.Text = 'Export Eigenvectors (.csv)';

            % Create SettingMenu
            app.SettingMenu = uimenu(app.UIFigure);
            app.SettingMenu.Text = 'Setting';

            % Create ThemeMenu
            app.ThemeMenu = uimenu(app.SettingMenu);
            app.ThemeMenu.Text = 'Theme';

            % Create WhiteMenu
            app.WhiteMenu = uimenu(app.ThemeMenu);
            app.WhiteMenu.MenuSelectedFcn = createCallbackFcn(app, @WhiteMenuSelected, true);
            app.WhiteMenu.Text = 'White';

            % Create BlueMenu
            app.BlueMenu = uimenu(app.ThemeMenu);
            app.BlueMenu.MenuSelectedFcn = createCallbackFcn(app, @BlueMenuSelected, true);
            app.BlueMenu.Text = 'Blue';

            % Create DarkMenu
            app.DarkMenu = uimenu(app.ThemeMenu);
            app.DarkMenu.MenuSelectedFcn = createCallbackFcn(app, @DarkMenuSelected, true);
            app.DarkMenu.Text = 'Dark';

            % Create DefaultMenu
            app.DefaultMenu = uimenu(app.ThemeMenu);
            app.DefaultMenu.MenuSelectedFcn = createCallbackFcn(app, @DefaultMenuSelected, true);
            app.DefaultMenu.Text = 'Default';

            % Create OwlsMenu
            app.OwlsMenu = uimenu(app.ThemeMenu);
            app.OwlsMenu.Text = 'Owls';

            % Create ShowmemoreMenu
            app.ShowmemoreMenu = uimenu(app.OwlsMenu);
            app.ShowmemoreMenu.MenuSelectedFcn = createCallbackFcn(app, @ShowmemoreMenuSelected, true);
            app.ShowmemoreMenu.Text = 'Show me more';

            % Create FlythemoffMenu
            app.FlythemoffMenu = uimenu(app.OwlsMenu);
            app.FlythemoffMenu.MenuSelectedFcn = createCallbackFcn(app, @FlythemoffMenuSelected, true);
            app.FlythemoffMenu.Text = 'Fly them off';

            % Create ResolutionMenu
            app.ResolutionMenu = uimenu(app.SettingMenu);
            app.ResolutionMenu.Text = 'Resolution';

            % Create SettoHighMenu
            app.SettoHighMenu = uimenu(app.ResolutionMenu);
            app.SettoHighMenu.MenuSelectedFcn = createCallbackFcn(app, @SettoHighMenuSelected, true);
            app.SettoHighMenu.Text = 'Set to High';

            % Create SettoMediumMenu
            app.SettoMediumMenu = uimenu(app.ResolutionMenu);
            app.SettoMediumMenu.MenuSelectedFcn = createCallbackFcn(app, @SettoMediumMenuSelected, true);
            app.SettoMediumMenu.Text = 'Set to Medium';

            % Create SettoLowMenu
            app.SettoLowMenu = uimenu(app.ResolutionMenu);
            app.SettoLowMenu.MenuSelectedFcn = createCallbackFcn(app, @SettoLowMenuSelected, true);
            app.SettoLowMenu.Text = 'Set to Low';

            % Create InfoonResolutionMenu
            app.InfoonResolutionMenu = uimenu(app.ResolutionMenu);
            app.InfoonResolutionMenu.MenuSelectedFcn = createCallbackFcn(app, @InfoonResolutionMenuSelected, true);
            app.InfoonResolutionMenu.Text = 'Info. on Resolution';

            % Create HiddenPotentialSettingsMenu
            app.HiddenPotentialSettingsMenu = uimenu(app.SettingMenu);
            app.HiddenPotentialSettingsMenu.Text = 'Hidden Potential Settings';

            % Create CoulmbicMenu
            app.CoulmbicMenu = uimenu(app.HiddenPotentialSettingsMenu);
            app.CoulmbicMenu.Text = 'Coulmbic ';

            % Create SettingMenu_4
            app.SettingMenu_4 = uimenu(app.CoulmbicMenu);
            app.SettingMenu_4.MenuSelectedFcn = createCallbackFcn(app, @SettingMenu_4Selected, true);
            app.SettingMenu_4.Text = 'Setting';

            % Create HelpMenu_4
            app.HelpMenu_4 = uimenu(app.CoulmbicMenu);
            app.HelpMenu_4.MenuSelectedFcn = createCallbackFcn(app, @HelpMenu_4Selected, true);
            app.HelpMenu_4.Text = 'Help';

            % Create HarmonicOscilatorMenu
            app.HarmonicOscilatorMenu = uimenu(app.HiddenPotentialSettingsMenu);
            app.HarmonicOscilatorMenu.Text = 'Harmonic Oscilator';

            % Create SettingMenu_3
            app.SettingMenu_3 = uimenu(app.HarmonicOscilatorMenu);
            app.SettingMenu_3.MenuSelectedFcn = createCallbackFcn(app, @SettingMenu_3Selected, true);
            app.SettingMenu_3.Text = 'Setting';

            % Create HelpMenu_3
            app.HelpMenu_3 = uimenu(app.HarmonicOscilatorMenu);
            app.HelpMenu_3.MenuSelectedFcn = createCallbackFcn(app, @HelpMenu_3Selected, true);
            app.HelpMenu_3.Text = 'Help';

            % Create MorseOscillatorMenu
            app.MorseOscillatorMenu = uimenu(app.HiddenPotentialSettingsMenu);
            app.MorseOscillatorMenu.MenuSelectedFcn = createCallbackFcn(app, @MorseOscillatorMenuSelected, true);
            app.MorseOscillatorMenu.Text = 'Morse Oscillator';

            % Create SettingMenu_2
            app.SettingMenu_2 = uimenu(app.MorseOscillatorMenu);
            app.SettingMenu_2.MenuSelectedFcn = createCallbackFcn(app, @SettingMenu_2Selected, true);
            app.SettingMenu_2.Text = 'Setting';

            % Create HelpMenu_2
            app.HelpMenu_2 = uimenu(app.MorseOscillatorMenu);
            app.HelpMenu_2.MenuSelectedFcn = createCallbackFcn(app, @HelpMenu_2Selected, true);
            app.HelpMenu_2.Text = 'Help';

            % Create AllowNegativeValuesofPotentialEnergyMenu
            app.AllowNegativeValuesofPotentialEnergyMenu = uimenu(app.HiddenPotentialSettingsMenu);
            app.AllowNegativeValuesofPotentialEnergyMenu.MenuSelectedFcn = createCallbackFcn(app, @AllowNegativeValuesofPotentialEnergyMenuSelected, true);
            app.AllowNegativeValuesofPotentialEnergyMenu.Text = 'Allow Negative Values of Potential Energy';

            % Create ShiftWavefunctionbyEnergyMenu
            app.ShiftWavefunctionbyEnergyMenu = uimenu(app.SettingMenu);
            app.ShiftWavefunctionbyEnergyMenu.MenuSelectedFcn = createCallbackFcn(app, @ShiftWavefunctionbyEnergyMenuSelected, true);
            app.ShiftWavefunctionbyEnergyMenu.Text = 'Shift Wavefunction by Energy';

            % Create EnableorDisableMenu
            app.EnableorDisableMenu = uimenu(app.ShiftWavefunctionbyEnergyMenu);
            app.EnableorDisableMenu.MenuSelectedFcn = createCallbackFcn(app, @EnableorDisableMenuSelected, true);
            app.EnableorDisableMenu.Text = 'Enable or Disable';

            % Create ShowEnergySpectrumMenu
            app.ShowEnergySpectrumMenu = uimenu(app.SettingMenu);
            app.ShowEnergySpectrumMenu.Text = 'Show Energy Spectrum';

            % Create SelectandPlotEigenvaluesMenu
            app.SelectandPlotEigenvaluesMenu = uimenu(app.ShowEnergySpectrumMenu);
            app.SelectandPlotEigenvaluesMenu.MenuSelectedFcn = createCallbackFcn(app, @SelectandPlotEigenvaluesMenuSelected, true);
            app.SelectandPlotEigenvaluesMenu.Text = 'Select and Plot Eigenvalues';

            % Create EnergySpectrumPlottingOptionsMenu
            app.EnergySpectrumPlottingOptionsMenu = uimenu(app.ShowEnergySpectrumMenu);
            app.EnergySpectrumPlottingOptionsMenu.Text = 'Energy Spectrum Plotting Options';

            % Create ChangePlotTypeMenu
            app.ChangePlotTypeMenu = uimenu(app.EnergySpectrumPlottingOptionsMenu);
            app.ChangePlotTypeMenu.MenuSelectedFcn = createCallbackFcn(app, @ChangePlotTypeMenuSelected, true);
            app.ChangePlotTypeMenu.Text = 'Change Plot Type';

            % Create NormalizationConstantMenu
            app.NormalizationConstantMenu = uimenu(app.SettingMenu);
            app.NormalizationConstantMenu.Text = 'Normalization Constant';

            % Create ChangeNormalizationConstantMenu
            app.ChangeNormalizationConstantMenu = uimenu(app.NormalizationConstantMenu);
            app.ChangeNormalizationConstantMenu.MenuSelectedFcn = createCallbackFcn(app, @ChangeNormalizationConstantMenuSelected, true);
            app.ChangeNormalizationConstantMenu.Text = 'Change Normalization Constant';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.UIFigure);
            app.HelpMenu.Text = 'Help';

            % Create DocumentationMenu
            app.DocumentationMenu = uimenu(app.HelpMenu);
            app.DocumentationMenu.MenuSelectedFcn = createCallbackFcn(app, @DocumentationMenuSelected, true);
            app.DocumentationMenu.Text = 'Documentation';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.HelpMenu);
            app.AboutMenu.MenuSelectedFcn = createCallbackFcn(app, @AboutMenuSelected, true);
            app.AboutMenu.Text = 'About';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.FontWeight = 'bold';
            app.UIAxes.Position = [3 2 669 616];

            % Create RunButton
            app.RunButton = uibutton(app.UIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Icon = 'OIP.erFIaMIASGYeeSlAESZltgAAAA.jpg';
            app.RunButton.BackgroundColor = [1 1 1];
            app.RunButton.FontSize = 20;
            app.RunButton.FontWeight = 'bold';
            app.RunButton.FontColor = [0.0078 0.0078 0.4196];
            app.RunButton.Position = [918 33 126 88];
            app.RunButton.Text = 'Run';

            % Create PeriodicBoundaryConditionDropDownLabel
            app.PeriodicBoundaryConditionDropDownLabel = uilabel(app.UIFigure);
            app.PeriodicBoundaryConditionDropDownLabel.BackgroundColor = [0.902 0.902 0.902];
            app.PeriodicBoundaryConditionDropDownLabel.HorizontalAlignment = 'center';
            app.PeriodicBoundaryConditionDropDownLabel.FontWeight = 'bold';
            app.PeriodicBoundaryConditionDropDownLabel.Position = [708 55 172 22];
            app.PeriodicBoundaryConditionDropDownLabel.Text = 'Periodic Boundary Condition';

            % Create PeriodicBoundaryConditionDropDown
            app.PeriodicBoundaryConditionDropDown = uidropdown(app.UIFigure);
            app.PeriodicBoundaryConditionDropDown.Items = {'Yes', 'No'};
            app.PeriodicBoundaryConditionDropDown.FontWeight = 'bold';
            app.PeriodicBoundaryConditionDropDown.BackgroundColor = [0.902 0.902 0.902];
            app.PeriodicBoundaryConditionDropDown.Position = [782 33 83 22];
            app.PeriodicBoundaryConditionDropDown.Value = 'No';

            % Create PotentialEnergyPanel
            app.PotentialEnergyPanel = uipanel(app.UIFigure);
            app.PotentialEnergyPanel.ForegroundColor = [0.4941 0.1843 0.5569];
            app.PotentialEnergyPanel.TitlePosition = 'centertop';
            app.PotentialEnergyPanel.Title = 'Potential Energy';
            app.PotentialEnergyPanel.FontWeight = 'bold';
            app.PotentialEnergyPanel.Position = [901 407 161 258];

            % Create FiniteSquareWellCheckBox
            app.FiniteSquareWellCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.FiniteSquareWellCheckBox.ValueChangedFcn = createCallbackFcn(app, @FiniteSquareWellCheckBoxValueChanged, true);
            app.FiniteSquareWellCheckBox.Text = 'Finite Square Well';
            app.FiniteSquareWellCheckBox.FontWeight = 'bold';
            app.FiniteSquareWellCheckBox.Position = [9 211 126 22];

            % Create FiniteSquareBarrierCheckBox
            app.FiniteSquareBarrierCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.FiniteSquareBarrierCheckBox.ValueChangedFcn = createCallbackFcn(app, @FiniteSquareBarrierCheckBoxValueChanged, true);
            app.FiniteSquareBarrierCheckBox.Text = 'Finite Square Barrier';
            app.FiniteSquareBarrierCheckBox.FontWeight = 'bold';
            app.FiniteSquareBarrierCheckBox.Position = [9 186 141 22];

            % Create CoulombPotentialCheckBox
            app.CoulombPotentialCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.CoulombPotentialCheckBox.Text = 'Coulomb Potential';
            app.CoulombPotentialCheckBox.FontWeight = 'bold';
            app.CoulombPotentialCheckBox.Position = [8 160 128 22];

            % Create StepPotentialCheckBox
            app.StepPotentialCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.StepPotentialCheckBox.ValueChangedFcn = createCallbackFcn(app, @StepPotentialCheckBoxValueChanged, true);
            app.StepPotentialCheckBox.Text = 'Step Potential';
            app.StepPotentialCheckBox.FontWeight = 'bold';
            app.StepPotentialCheckBox.Position = [8 134 102 22];

            % Create HarmonicOscillatorCheckBox
            app.HarmonicOscillatorCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.HarmonicOscillatorCheckBox.ValueChangedFcn = createCallbackFcn(app, @HarmonicOscillatorCheckBoxValueChanged, true);
            app.HarmonicOscillatorCheckBox.Text = 'Harmonic Oscillator';
            app.HarmonicOscillatorCheckBox.FontWeight = 'bold';
            app.HarmonicOscillatorCheckBox.Position = [8 109 136 22];

            % Create MorseOscillatorCheckBox
            app.MorseOscillatorCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.MorseOscillatorCheckBox.ValueChangedFcn = createCallbackFcn(app, @MorseOscillatorCheckBoxValueChanged, true);
            app.MorseOscillatorCheckBox.Text = 'Morse Oscillator';
            app.MorseOscillatorCheckBox.FontWeight = 'bold';
            app.MorseOscillatorCheckBox.Position = [8 83 116 22];

            % Create RandomPotentialCheckBox
            app.RandomPotentialCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.RandomPotentialCheckBox.ValueChangedFcn = createCallbackFcn(app, @RandomPotentialCheckBoxValueChanged, true);
            app.RandomPotentialCheckBox.Text = 'Random Potential';
            app.RandomPotentialCheckBox.FontWeight = 'bold';
            app.RandomPotentialCheckBox.Position = [8 31 124 22];

            % Create DoubleWellCheckBox
            app.DoubleWellCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.DoubleWellCheckBox.ValueChangedFcn = createCallbackFcn(app, @DoubleWellCheckBoxValueChanged, true);
            app.DoubleWellCheckBox.Text = 'Double Well';
            app.DoubleWellCheckBox.FontWeight = 'bold';
            app.DoubleWellCheckBox.Position = [8 57 90 22];

            % Create ImportedPotentialCheckBox
            app.ImportedPotentialCheckBox = uicheckbox(app.PotentialEnergyPanel);
            app.ImportedPotentialCheckBox.ValueChangedFcn = createCallbackFcn(app, @ImportedPotentialCheckBoxValueChanged, true);
            app.ImportedPotentialCheckBox.Text = 'Imported Potential';
            app.ImportedPotentialCheckBox.FontWeight = 'bold';
            app.ImportedPotentialCheckBox.Position = [8 5 127 22];

            % Create EigenstateVisualizationPanel
            app.EigenstateVisualizationPanel = uipanel(app.UIFigure);
            app.EigenstateVisualizationPanel.TitlePosition = 'centertop';
            app.EigenstateVisualizationPanel.Title = 'Eigenstate Visualization';
            app.EigenstateVisualizationPanel.FontWeight = 'bold';
            app.EigenstateVisualizationPanel.Position = [900 138 164 212];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.EigenstateVisualizationPanel);
            app.TabGroup.Position = [1 2 161 188];

            % Create SingleTab
            app.SingleTab = uitab(app.TabGroup);
            app.SingleTab.Title = 'Single';

            % Create EigenStateQuantumNumberLabel
            app.EigenStateQuantumNumberLabel = uilabel(app.SingleTab);
            app.EigenStateQuantumNumberLabel.HorizontalAlignment = 'center';
            app.EigenStateQuantumNumberLabel.FontWeight = 'bold';
            app.EigenStateQuantumNumberLabel.Position = [27 123 107 28];
            app.EigenStateQuantumNumberLabel.Text = {'EigenState '; 'Quantum Number'};

            % Create EigenStateQuantumNumberSpinner
            app.EigenStateQuantumNumberSpinner = uispinner(app.SingleTab);
            app.EigenStateQuantumNumberSpinner.Limits = [1 1000];
            app.EigenStateQuantumNumberSpinner.HorizontalAlignment = 'center';
            app.EigenStateQuantumNumberSpinner.FontWeight = 'bold';
            app.EigenStateQuantumNumberSpinner.Position = [17 85 128 22];
            app.EigenStateQuantumNumberSpinner.Value = 1;

            % Create Wavefunction2Label
            app.Wavefunction2Label = uilabel(app.SingleTab);
            app.Wavefunction2Label.HorizontalAlignment = 'center';
            app.Wavefunction2Label.FontWeight = 'bold';
            app.Wavefunction2Label.FontColor = [0.0314 0.0039 0.2706];
            app.Wavefunction2Label.Position = [30 56 105 22];
            app.Wavefunction2Label.Text = '(Wavefunction)^2';

            % Create Wavefunction2Switch
            app.Wavefunction2Switch = uiswitch(app.SingleTab, 'slider');
            app.Wavefunction2Switch.ValueChangedFcn = createCallbackFcn(app, @Wavefunction2SwitchValueChanged, true);
            app.Wavefunction2Switch.FontWeight = 'bold';
            app.Wavefunction2Switch.FontColor = [0.0314 0.0039 0.2706];
            app.Wavefunction2Switch.Position = [56 30 54 24];

            % Create MultiStateVisualizationmodeLabel
            app.MultiStateVisualizationmodeLabel = uilabel(app.SingleTab);
            app.MultiStateVisualizationmodeLabel.HorizontalAlignment = 'center';
            app.MultiStateVisualizationmodeLabel.FontWeight = 'bold';
            app.MultiStateVisualizationmodeLabel.Visible = 'off';
            app.MultiStateVisualizationmodeLabel.Position = [27 13 114 28];
            app.MultiStateVisualizationmodeLabel.Text = {'Multi-State '; 'Visualization mode'};

            % Create MultiTab
            app.MultiTab = uitab(app.TabGroup);
            app.MultiTab.Title = 'Multi';

            % Create MultiSwitch
            app.MultiSwitch = uiswitch(app.MultiTab, 'slider');
            app.MultiSwitch.ValueChangedFcn = createCallbackFcn(app, @MultiSwitchValueChanged, true);
            app.MultiSwitch.FontWeight = 'bold';
            app.MultiSwitch.Position = [64 130 45 20];

            % Create StateSpinnerLabel
            app.StateSpinnerLabel = uilabel(app.MultiTab);
            app.StateSpinnerLabel.HorizontalAlignment = 'right';
            app.StateSpinnerLabel.FontWeight = 'bold';
            app.StateSpinnerLabel.Position = [-1 100 45 22];
            app.StateSpinnerLabel.Text = 'State #';

            % Create StateSpinner
            app.StateSpinner = uispinner(app.MultiTab);
            app.StateSpinner.Limits = [0 1000];
            app.StateSpinner.FontWeight = 'bold';
            app.StateSpinner.Position = [59 100 100 22];
            app.StateSpinner.Value = 1;

            % Create StateSpinner_2Label
            app.StateSpinner_2Label = uilabel(app.MultiTab);
            app.StateSpinner_2Label.HorizontalAlignment = 'right';
            app.StateSpinner_2Label.FontWeight = 'bold';
            app.StateSpinner_2Label.FontColor = [0.851 0.3255 0.098];
            app.StateSpinner_2Label.Position = [-1 69 45 22];
            app.StateSpinner_2Label.Text = 'State #';

            % Create StateSpinner_2
            app.StateSpinner_2 = uispinner(app.MultiTab);
            app.StateSpinner_2.Limits = [0 1000];
            app.StateSpinner_2.FontWeight = 'bold';
            app.StateSpinner_2.FontColor = [0.851 0.3255 0.098];
            app.StateSpinner_2.Position = [59 69 100 22];
            app.StateSpinner_2.Value = 2;

            % Create StateSpinner_3Label
            app.StateSpinner_3Label = uilabel(app.MultiTab);
            app.StateSpinner_3Label.HorizontalAlignment = 'right';
            app.StateSpinner_3Label.FontWeight = 'bold';
            app.StateSpinner_3Label.FontColor = [0 0 1];
            app.StateSpinner_3Label.Position = [-1 38 45 22];
            app.StateSpinner_3Label.Text = 'State #';

            % Create StateSpinner_3
            app.StateSpinner_3 = uispinner(app.MultiTab);
            app.StateSpinner_3.Limits = [0 1000];
            app.StateSpinner_3.FontWeight = 'bold';
            app.StateSpinner_3.FontColor = [0 0 1];
            app.StateSpinner_3.Position = [59 38 100 22];

            % Create StateSpinner_4Label
            app.StateSpinner_4Label = uilabel(app.MultiTab);
            app.StateSpinner_4Label.HorizontalAlignment = 'right';
            app.StateSpinner_4Label.FontWeight = 'bold';
            app.StateSpinner_4Label.FontColor = [0.3922 0.8314 0.0745];
            app.StateSpinner_4Label.Position = [-1 5 45 22];
            app.StateSpinner_4Label.Text = 'State #';

            % Create StateSpinner_4
            app.StateSpinner_4 = uispinner(app.MultiTab);
            app.StateSpinner_4.Limits = [0 1000];
            app.StateSpinner_4.FontWeight = 'bold';
            app.StateSpinner_4.FontColor = [0.3922 0.8314 0.0745];
            app.StateSpinner_4.Position = [59 5 100 22];

            % Create Image
            app.Image = uiimage(app.UIFigure);
            app.Image.ImageClickedFcn = createCallbackFcn(app, @ImageClicked, true);
            app.Image.Position = [3 616 428 61];
            app.Image.ImageSource = 'Rice_Horizontal_Blue_1.png';

            % Create MassDropDownLabel
            app.MassDropDownLabel = uilabel(app.UIFigure);
            app.MassDropDownLabel.BackgroundColor = [0.902 0.902 0.902];
            app.MassDropDownLabel.FontWeight = 'bold';
            app.MassDropDownLabel.FontColor = [0 0 1];
            app.MassDropDownLabel.Position = [903 373 40 22];
            app.MassDropDownLabel.Text = 'Mass:';

            % Create MassDropDown
            app.MassDropDown = uidropdown(app.UIFigure);
            app.MassDropDown.Items = {'Electron', 'Proton', 'Helium', 'H2'};
            app.MassDropDown.FontWeight = 'bold';
            app.MassDropDown.FontColor = [0 0 1];
            app.MassDropDown.BackgroundColor = [0.902 0.902 0.902];
            app.MassDropDown.Position = [947 373 110 22];
            app.MassDropDown.Value = 'Electron';

            % Create PotentialModificationPanel
            app.PotentialModificationPanel = uipanel(app.UIFigure);
            app.PotentialModificationPanel.TitlePosition = 'centertop';
            app.PotentialModificationPanel.Title = 'Potential Modification';
            app.PotentialModificationPanel.BackgroundColor = [0.902 0.902 0.902];
            app.PotentialModificationPanel.FontWeight = 'bold';
            app.PotentialModificationPanel.Position = [679 88 218 577];

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.PotentialModificationPanel);
            app.TabGroup2.Position = [4 397 213 156];

            % Create WellTab
            app.WellTab = uitab(app.TabGroup2);
            app.WellTab.Title = 'Well';

            % Create WellDepthEditFieldLabel
            app.WellDepthEditFieldLabel = uilabel(app.WellTab);
            app.WellDepthEditFieldLabel.HorizontalAlignment = 'center';
            app.WellDepthEditFieldLabel.FontWeight = 'bold';
            app.WellDepthEditFieldLabel.Position = [46 38 67 22];
            app.WellDepthEditFieldLabel.Text = 'Well Depth';

            % Create WellDepthEditField
            app.WellDepthEditField = uieditfield(app.WellTab, 'numeric');
            app.WellDepthEditField.Limits = [0 Inf];
            app.WellDepthEditField.HorizontalAlignment = 'center';
            app.WellDepthEditField.FontWeight = 'bold';
            app.WellDepthEditField.Position = [126 38 68 22];
            app.WellDepthEditField.Value = 2;

            % Create LeftWallPositionEditFieldLabel
            app.LeftWallPositionEditFieldLabel = uilabel(app.WellTab);
            app.LeftWallPositionEditFieldLabel.HorizontalAlignment = 'center';
            app.LeftWallPositionEditFieldLabel.FontWeight = 'bold';
            app.LeftWallPositionEditFieldLabel.Position = [12 98 106 22];
            app.LeftWallPositionEditFieldLabel.Text = 'Left Wall Position';

            % Create LeftWallPositionEditField
            app.LeftWallPositionEditField = uieditfield(app.WellTab, 'numeric');
            app.LeftWallPositionEditField.Limits = [0 1];
            app.LeftWallPositionEditField.ValueChangedFcn = createCallbackFcn(app, @LeftWallPositionEditFieldValueChanged, true);
            app.LeftWallPositionEditField.HorizontalAlignment = 'center';
            app.LeftWallPositionEditField.FontWeight = 'bold';
            app.LeftWallPositionEditField.Position = [126 98 68 22];
            app.LeftWallPositionEditField.Value = 0.4;

            % Create RightWallPositionEditFieldLabel
            app.RightWallPositionEditFieldLabel = uilabel(app.WellTab);
            app.RightWallPositionEditFieldLabel.HorizontalAlignment = 'center';
            app.RightWallPositionEditFieldLabel.FontWeight = 'bold';
            app.RightWallPositionEditFieldLabel.Position = [9 68 114 22];
            app.RightWallPositionEditFieldLabel.Text = 'Right Wall Position';

            % Create RightWallPositionEditField
            app.RightWallPositionEditField = uieditfield(app.WellTab, 'numeric');
            app.RightWallPositionEditField.Limits = [0 1];
            app.RightWallPositionEditField.ValueChangedFcn = createCallbackFcn(app, @RightWallPositionEditFieldValueChanged, true);
            app.RightWallPositionEditField.HorizontalAlignment = 'center';
            app.RightWallPositionEditField.FontWeight = 'bold';
            app.RightWallPositionEditField.Position = [126 68 68 22];
            app.RightWallPositionEditField.Value = 0.6;

            % Create ResettoDefaultButton
            app.ResettoDefaultButton = uibutton(app.WellTab, 'push');
            app.ResettoDefaultButton.ButtonPushedFcn = createCallbackFcn(app, @ResettoDefaultButtonPushed, true);
            app.ResettoDefaultButton.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton.FontWeight = 'bold';
            app.ResettoDefaultButton.Position = [53 8 107 22];
            app.ResettoDefaultButton.Text = 'Reset to Default';

            % Create BarrierTab
            app.BarrierTab = uitab(app.TabGroup2);
            app.BarrierTab.Title = 'Barrier';

            % Create ResettoDefaultButton_2
            app.ResettoDefaultButton_2 = uibutton(app.BarrierTab, 'push');
            app.ResettoDefaultButton_2.ButtonPushedFcn = createCallbackFcn(app, @ResettoDefaultButton_2Pushed, true);
            app.ResettoDefaultButton_2.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton_2.FontWeight = 'bold';
            app.ResettoDefaultButton_2.Position = [53 10 107 22];
            app.ResettoDefaultButton_2.Text = 'Reset to Default';

            % Create LeftWallPositionEditField_2Label
            app.LeftWallPositionEditField_2Label = uilabel(app.BarrierTab);
            app.LeftWallPositionEditField_2Label.HorizontalAlignment = 'right';
            app.LeftWallPositionEditField_2Label.FontWeight = 'bold';
            app.LeftWallPositionEditField_2Label.Position = [12 98 106 22];
            app.LeftWallPositionEditField_2Label.Text = 'Left Wall Position';

            % Create LeftWallPositionEditField_2
            app.LeftWallPositionEditField_2 = uieditfield(app.BarrierTab, 'numeric');
            app.LeftWallPositionEditField_2.Limits = [0 1];
            app.LeftWallPositionEditField_2.ValueChangedFcn = createCallbackFcn(app, @LeftWallPositionEditField_2ValueChanged, true);
            app.LeftWallPositionEditField_2.HorizontalAlignment = 'center';
            app.LeftWallPositionEditField_2.FontWeight = 'bold';
            app.LeftWallPositionEditField_2.Position = [132 98 67 22];
            app.LeftWallPositionEditField_2.Value = 0.45;

            % Create RightWallPositionEditField_2Label
            app.RightWallPositionEditField_2Label = uilabel(app.BarrierTab);
            app.RightWallPositionEditField_2Label.HorizontalAlignment = 'right';
            app.RightWallPositionEditField_2Label.FontWeight = 'bold';
            app.RightWallPositionEditField_2Label.Position = [5 68 114 22];
            app.RightWallPositionEditField_2Label.Text = 'Right Wall Position';

            % Create RightWallPositionEditField_2
            app.RightWallPositionEditField_2 = uieditfield(app.BarrierTab, 'numeric');
            app.RightWallPositionEditField_2.Limits = [0 1];
            app.RightWallPositionEditField_2.ValueChangedFcn = createCallbackFcn(app, @RightWallPositionEditField_2ValueChanged, true);
            app.RightWallPositionEditField_2.HorizontalAlignment = 'center';
            app.RightWallPositionEditField_2.FontWeight = 'bold';
            app.RightWallPositionEditField_2.Position = [132 68 67 22];
            app.RightWallPositionEditField_2.Value = 0.55;

            % Create BarrierHeightEditFieldLabel
            app.BarrierHeightEditFieldLabel = uilabel(app.BarrierTab);
            app.BarrierHeightEditFieldLabel.HorizontalAlignment = 'right';
            app.BarrierHeightEditFieldLabel.FontWeight = 'bold';
            app.BarrierHeightEditFieldLabel.Position = [31 37 86 22];
            app.BarrierHeightEditFieldLabel.Text = 'Barrier Height';

            % Create BarrierHeightEditField
            app.BarrierHeightEditField = uieditfield(app.BarrierTab, 'numeric');
            app.BarrierHeightEditField.Limits = [0 Inf];
            app.BarrierHeightEditField.HorizontalAlignment = 'center';
            app.BarrierHeightEditField.FontWeight = 'bold';
            app.BarrierHeightEditField.Position = [132 37 67 22];
            app.BarrierHeightEditField.Value = 1;

            % Create StepTab
            app.StepTab = uitab(app.TabGroup2);
            app.StepTab.Title = 'Step';

            % Create WallPositionEditFieldLabel
            app.WallPositionEditFieldLabel = uilabel(app.StepTab);
            app.WallPositionEditFieldLabel.HorizontalAlignment = 'right';
            app.WallPositionEditFieldLabel.FontWeight = 'bold';
            app.WallPositionEditFieldLabel.Position = [21 39 80 22];
            app.WallPositionEditFieldLabel.Text = 'Wall Position';

            % Create WallPositionEditField
            app.WallPositionEditField = uieditfield(app.StepTab, 'numeric');
            app.WallPositionEditField.Limits = [0 1];
            app.WallPositionEditField.HorizontalAlignment = 'center';
            app.WallPositionEditField.FontWeight = 'bold';
            app.WallPositionEditField.Position = [120 39 67 22];
            app.WallPositionEditField.Value = 0.5;

            % Create RightPotentialEditFieldLabel
            app.RightPotentialEditFieldLabel = uilabel(app.StepTab);
            app.RightPotentialEditFieldLabel.HorizontalAlignment = 'right';
            app.RightPotentialEditFieldLabel.FontWeight = 'bold';
            app.RightPotentialEditFieldLabel.Position = [16 95 90 22];
            app.RightPotentialEditFieldLabel.Text = 'Right Potential';

            % Create RightPotentialEditField
            app.RightPotentialEditField = uieditfield(app.StepTab, 'numeric');
            app.RightPotentialEditField.Limits = [0 Inf];
            app.RightPotentialEditField.ValueChangedFcn = createCallbackFcn(app, @RightPotentialEditFieldValueChanged, true);
            app.RightPotentialEditField.HorizontalAlignment = 'center';
            app.RightPotentialEditField.FontWeight = 'bold';
            app.RightPotentialEditField.Position = [120 95 67 22];
            app.RightPotentialEditField.Value = 1;

            % Create LeftPotentialEditFieldLabel
            app.LeftPotentialEditFieldLabel = uilabel(app.StepTab);
            app.LeftPotentialEditFieldLabel.HorizontalAlignment = 'right';
            app.LeftPotentialEditFieldLabel.FontWeight = 'bold';
            app.LeftPotentialEditFieldLabel.Position = [21 67 82 22];
            app.LeftPotentialEditFieldLabel.Text = 'Left Potential';

            % Create LeftPotentialEditField
            app.LeftPotentialEditField = uieditfield(app.StepTab, 'numeric');
            app.LeftPotentialEditField.Limits = [0 Inf];
            app.LeftPotentialEditField.HorizontalAlignment = 'center';
            app.LeftPotentialEditField.FontWeight = 'bold';
            app.LeftPotentialEditField.Position = [120 67 67 22];
            app.LeftPotentialEditField.Value = 0.1;

            % Create ResettoDefaultButton_3
            app.ResettoDefaultButton_3 = uibutton(app.StepTab, 'push');
            app.ResettoDefaultButton_3.ButtonPushedFcn = createCallbackFcn(app, @ResettoDefaultButton_3Pushed, true);
            app.ResettoDefaultButton_3.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton_3.FontWeight = 'bold';
            app.ResettoDefaultButton_3.Position = [53 10 107 22];
            app.ResettoDefaultButton_3.Text = 'Reset to Default';

            % Create TabGroup3
            app.TabGroup3 = uitabgroup(app.PotentialModificationPanel);
            app.TabGroup3.Position = [2 239 213 150];

            % Create HarmonicTab
            app.HarmonicTab = uitab(app.TabGroup3);
            app.HarmonicTab.Title = 'Harmonic';

            % Create ResettoDefaultButton_4
            app.ResettoDefaultButton_4 = uibutton(app.HarmonicTab, 'push');
            app.ResettoDefaultButton_4.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton_4.FontWeight = 'bold';
            app.ResettoDefaultButton_4.Position = [55 25 107 22];
            app.ResettoDefaultButton_4.Text = 'Reset to Default';

            % Create EquPositionEditFieldLabel
            app.EquPositionEditFieldLabel = uilabel(app.HarmonicTab);
            app.EquPositionEditFieldLabel.HorizontalAlignment = 'right';
            app.EquPositionEditFieldLabel.FontWeight = 'bold';
            app.EquPositionEditFieldLabel.Position = [31 65 82 22];
            app.EquPositionEditFieldLabel.Text = 'Equ. Position';

            % Create EquPositionEditField
            app.EquPositionEditField = uieditfield(app.HarmonicTab, 'numeric');
            app.EquPositionEditField.Limits = [0 1];
            app.EquPositionEditField.HorizontalAlignment = 'center';
            app.EquPositionEditField.FontWeight = 'bold';
            app.EquPositionEditField.Position = [127 65 67 22];
            app.EquPositionEditField.Value = 0.5;

            % Create ForceConstantEditFieldLabel
            app.ForceConstantEditFieldLabel = uilabel(app.HarmonicTab);
            app.ForceConstantEditFieldLabel.HorizontalAlignment = 'right';
            app.ForceConstantEditFieldLabel.FontWeight = 'bold';
            app.ForceConstantEditFieldLabel.Position = [19 95 94 22];
            app.ForceConstantEditFieldLabel.Text = 'Force Constant';

            % Create ForceConstantEditField
            app.ForceConstantEditField = uieditfield(app.HarmonicTab, 'numeric');
            app.ForceConstantEditField.Limits = [0 Inf];
            app.ForceConstantEditField.HorizontalAlignment = 'center';
            app.ForceConstantEditField.FontWeight = 'bold';
            app.ForceConstantEditField.Position = [127 95 67 22];
            app.ForceConstantEditField.Value = 0.004;

            % Create MorseTab
            app.MorseTab = uitab(app.TabGroup3);
            app.MorseTab.Title = 'Morse';

            % Create ResettoDefaultButton_5
            app.ResettoDefaultButton_5 = uibutton(app.MorseTab, 'push');
            app.ResettoDefaultButton_5.ButtonPushedFcn = createCallbackFcn(app, @ResettoDefaultButton_5Pushed, true);
            app.ResettoDefaultButton_5.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton_5.FontWeight = 'bold';
            app.ResettoDefaultButton_5.Position = [55 10 107 22];
            app.ResettoDefaultButton_5.Text = 'Reset to Default';

            % Create ForceConstantEditField_2Label
            app.ForceConstantEditField_2Label = uilabel(app.MorseTab);
            app.ForceConstantEditField_2Label.HorizontalAlignment = 'right';
            app.ForceConstantEditField_2Label.FontWeight = 'bold';
            app.ForceConstantEditField_2Label.Position = [15 94 94 22];
            app.ForceConstantEditField_2Label.Text = 'Force Constant';

            % Create ForceConstantEditField_2
            app.ForceConstantEditField_2 = uieditfield(app.MorseTab, 'numeric');
            app.ForceConstantEditField_2.Limits = [0 Inf];
            app.ForceConstantEditField_2.HorizontalAlignment = 'center';
            app.ForceConstantEditField_2.FontWeight = 'bold';
            app.ForceConstantEditField_2.Position = [128 94 67 22];
            app.ForceConstantEditField_2.Value = 0.363;

            % Create EquPositionEditField_2Label
            app.EquPositionEditField_2Label = uilabel(app.MorseTab);
            app.EquPositionEditField_2Label.HorizontalAlignment = 'right';
            app.EquPositionEditField_2Label.FontWeight = 'bold';
            app.EquPositionEditField_2Label.Position = [27 39 82 22];
            app.EquPositionEditField_2Label.Text = 'Equ. Position';

            % Create EquPositionEditField_2
            app.EquPositionEditField_2 = uieditfield(app.MorseTab, 'numeric');
            app.EquPositionEditField_2.Limits = [0 1];
            app.EquPositionEditField_2.HorizontalAlignment = 'center';
            app.EquPositionEditField_2.FontWeight = 'bold';
            app.EquPositionEditField_2.Position = [128 39 67 22];
            app.EquPositionEditField_2.Value = 0.5;

            % Create DissociationEnergyEditFieldLabel
            app.DissociationEnergyEditFieldLabel = uilabel(app.MorseTab);
            app.DissociationEnergyEditFieldLabel.HorizontalAlignment = 'right';
            app.DissociationEnergyEditFieldLabel.FontWeight = 'bold';
            app.DissociationEnergyEditFieldLabel.Position = [3 67 121 22];
            app.DissociationEnergyEditFieldLabel.Text = 'Dissociation Energy';

            % Create DissociationEnergyEditField
            app.DissociationEnergyEditField = uieditfield(app.MorseTab, 'numeric');
            app.DissociationEnergyEditField.Limits = [0 Inf];
            app.DissociationEnergyEditField.HorizontalAlignment = 'center';
            app.DissociationEnergyEditField.FontWeight = 'bold';
            app.DissociationEnergyEditField.Position = [128 67 67 22];
            app.DissociationEnergyEditField.Value = 0.1816;

            % Create RandomTab_2
            app.RandomTab_2 = uitab(app.TabGroup3);
            app.RandomTab_2.Title = 'Random';

            % Create NoiseWeightEditFieldLabel
            app.NoiseWeightEditFieldLabel = uilabel(app.RandomTab_2);
            app.NoiseWeightEditFieldLabel.HorizontalAlignment = 'center';
            app.NoiseWeightEditFieldLabel.FontWeight = 'bold';
            app.NoiseWeightEditFieldLabel.Position = [64 85 81 22];
            app.NoiseWeightEditFieldLabel.Text = 'Noise Weight';

            % Create NoiseWeightEditField
            app.NoiseWeightEditField = uieditfield(app.RandomTab_2, 'numeric');
            app.NoiseWeightEditField.Limits = [0 Inf];
            app.NoiseWeightEditField.HorizontalAlignment = 'center';
            app.NoiseWeightEditField.FontWeight = 'bold';
            app.NoiseWeightEditField.Position = [57 59 100 22];
            app.NoiseWeightEditField.Value = 0.1;

            % Create ResettoDefaultButton_7
            app.ResettoDefaultButton_7 = uibutton(app.RandomTab_2, 'push');
            app.ResettoDefaultButton_7.ButtonPushedFcn = createCallbackFcn(app, @ResettoDefaultButton_7Pushed, true);
            app.ResettoDefaultButton_7.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton_7.FontWeight = 'bold';
            app.ResettoDefaultButton_7.Position = [55 14 107 22];
            app.ResettoDefaultButton_7.Text = 'Reset to Default';

            % Create TabGroup4
            app.TabGroup4 = uitabgroup(app.PotentialModificationPanel);
            app.TabGroup4.Position = [5 9 209 222];

            % Create DoubelWellTab
            app.DoubelWellTab = uitab(app.TabGroup4);
            app.DoubelWellTab.Title = 'Doubel Well';

            % Create PEcaxx1xx2xx3xx4Label
            app.PEcaxx1xx2xx3xx4Label = uilabel(app.DoubelWellTab);
            app.PEcaxx1xx2xx3xx4Label.FontWeight = 'bold';
            app.PEcaxx1xx2xx3xx4Label.FontColor = [0.1608 0 0.3216];
            app.PEcaxx1xx2xx3xx4Label.Position = [16 166 187 22];
            app.PEcaxx1xx2xx3xx4Label.Text = 'PE=c.(a.x-x1).(x-x2).(x-x3).(x-x4)';

            % Create aEditFieldLabel
            app.aEditFieldLabel = uilabel(app.DoubelWellTab);
            app.aEditFieldLabel.HorizontalAlignment = 'center';
            app.aEditFieldLabel.FontWeight = 'bold';
            app.aEditFieldLabel.Position = [116 134 25 22];
            app.aEditFieldLabel.Text = 'a';

            % Create aEditField
            app.aEditField = uieditfield(app.DoubelWellTab, 'numeric');
            app.aEditField.Position = [150 136 50 22];
            app.aEditField.Value = 1;

            % Create x1EditFieldLabel
            app.x1EditFieldLabel = uilabel(app.DoubelWellTab);
            app.x1EditFieldLabel.HorizontalAlignment = 'center';
            app.x1EditFieldLabel.FontWeight = 'bold';
            app.x1EditFieldLabel.Position = [12 107 25 22];
            app.x1EditFieldLabel.Text = 'x1';

            % Create x1EditField
            app.x1EditField = uieditfield(app.DoubelWellTab, 'numeric');
            app.x1EditField.Limits = [0 1];
            app.x1EditField.Position = [48 107 50 22];
            app.x1EditField.Value = 0.02;

            % Create x2EditFieldLabel
            app.x2EditFieldLabel = uilabel(app.DoubelWellTab);
            app.x2EditFieldLabel.HorizontalAlignment = 'center';
            app.x2EditFieldLabel.FontWeight = 'bold';
            app.x2EditFieldLabel.Position = [116 107 25 22];
            app.x2EditFieldLabel.Text = 'x2';

            % Create x2EditField
            app.x2EditField = uieditfield(app.DoubelWellTab, 'numeric');
            app.x2EditField.Limits = [0 1];
            app.x2EditField.Position = [150 107 52 22];
            app.x2EditField.Value = 0.35;

            % Create x3EditFieldLabel
            app.x3EditFieldLabel = uilabel(app.DoubelWellTab);
            app.x3EditFieldLabel.HorizontalAlignment = 'center';
            app.x3EditFieldLabel.FontWeight = 'bold';
            app.x3EditFieldLabel.Position = [12 81 25 22];
            app.x3EditFieldLabel.Text = 'x3';

            % Create x3EditField
            app.x3EditField = uieditfield(app.DoubelWellTab, 'numeric');
            app.x3EditField.Limits = [0 1];
            app.x3EditField.Position = [48 81 50 22];
            app.x3EditField.Value = 0.65;

            % Create x4EditFieldLabel
            app.x4EditFieldLabel = uilabel(app.DoubelWellTab);
            app.x4EditFieldLabel.HorizontalAlignment = 'center';
            app.x4EditFieldLabel.FontWeight = 'bold';
            app.x4EditFieldLabel.Position = [116 82 25 22];
            app.x4EditFieldLabel.Text = 'x4';

            % Create x4EditField
            app.x4EditField = uieditfield(app.DoubelWellTab, 'numeric');
            app.x4EditField.Limits = [0 1];
            app.x4EditField.Position = [150 82 52 22];
            app.x4EditField.Value = 0.98;

            % Create ResettoDefaultButton_6
            app.ResettoDefaultButton_6 = uibutton(app.DoubelWellTab, 'push');
            app.ResettoDefaultButton_6.ButtonPushedFcn = createCallbackFcn(app, @ResettoDefaultButton_6Pushed, true);
            app.ResettoDefaultButton_6.BackgroundColor = [0.8 0.8 0.8];
            app.ResettoDefaultButton_6.FontWeight = 'bold';
            app.ResettoDefaultButton_6.Position = [52 50 107 22];
            app.ResettoDefaultButton_6.Text = 'Reset to Default';

            % Create cEditFieldLabel
            app.cEditFieldLabel = uilabel(app.DoubelWellTab);
            app.cEditFieldLabel.HorizontalAlignment = 'center';
            app.cEditFieldLabel.FontWeight = 'bold';
            app.cEditFieldLabel.Position = [12 137 25 22];
            app.cEditFieldLabel.Text = 'c';

            % Create cEditField
            app.cEditField = uieditfield(app.DoubelWellTab, 'numeric');
            app.cEditField.Enable = 'off';
            app.cEditField.Position = [48 136 50 22];
            app.cEditField.Value = 1;

            % Create ImportTab
            app.ImportTab = uitab(app.TabGroup4);
            app.ImportTab.Title = 'Import';

            % Create ImportPotentialButton
            app.ImportPotentialButton = uibutton(app.ImportTab, 'state');
            app.ImportPotentialButton.ValueChangedFcn = createCallbackFcn(app, @ImportPotentialButtonValueChanged, true);
            app.ImportPotentialButton.Text = 'Import Potential';
            app.ImportPotentialButton.BackgroundColor = [0.902 0.902 0.902];
            app.ImportPotentialButton.FontWeight = 'bold';
            app.ImportPotentialButton.Position = [45 106 107 22];

            % Create IMPRTPETextArea
            app.IMPRTPETextArea = uitextarea(app.ImportTab);
            app.IMPRTPETextArea.Editable = 'off';
            app.IMPRTPETextArea.Position = [11 134 185 55];

            % Create IMPRT_PELamp
            app.IMPRT_PELamp = uilamp(app.ImportTab);
            app.IMPRT_PELamp.Position = [11 107 20 20];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.ImportTab);
            title(app.UIAxes2, '')
            xlabel(app.UIAxes2, '')
            ylabel(app.UIAxes2, '')
            app.UIAxes2.FontSize = 1;
            app.UIAxes2.XTick = [];
            app.UIAxes2.XTickLabel = '';
            app.UIAxes2.YTick = [];
            app.UIAxes2.YTickLabel = '';
            app.UIAxes2.Position = [5 8 197 90];

            % Create Button
            app.Button = uibutton(app.ImportTab, 'state');
            app.Button.ValueChangedFcn = createCallbackFcn(app, @ButtonValueChanged, true);
            app.Button.Text = '?';
            app.Button.FontWeight = 'bold';
            app.Button.Position = [168 107 25 20];

            % Create Image2
            app.Image2 = uiimage(app.UIFigure);
            app.Image2.Position = [589 614 89 63];
            app.Image2.ImageSource = 'OIP.1A9GpaTOeUi6sjHi1BPEcQHaHl.jpg';

            % Create Image3
            app.Image3 = uiimage(app.UIFigure);
            app.Image3.Position = [429 614 61 63];
            app.Image3.ImageSource = 'OIP.Bn1I1NNo-QbkUD3n7KvGtgHaHm.jpg';

            % Create Image4
            app.Image4 = uiimage(app.UIFigure);
            app.Image4.Position = [498 614 100 63];
            app.Image4.ImageSource = 'R.ab66803f1bf1a9bf2bddca5209202fa1.png';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SchrodingerWave

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end