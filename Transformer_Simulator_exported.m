classdef Transformer_Simulator_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        GridLayout                     matlab.ui.container.GridLayout
        LeftPanel                      matlab.ui.container.Panel
        TransformerRatingSinKVAEditFieldLabel  matlab.ui.control.Label
        TransformerRatingSinKVAEditField  matlab.ui.control.NumericEditField
        PrimaryVoltageLVinVLabel       matlab.ui.control.Label
        PrimaryVoltageLVinVEditField   matlab.ui.control.NumericEditField
        SecondaryVoltageHVinKVLabel    matlab.ui.control.Label
        SecondaryVoltageHVinKVEditField  matlab.ui.control.NumericEditField
        CurrentDensityinAmm2EditFieldLabel  matlab.ui.control.Label
        CurrentDensityinAmm2EditField  matlab.ui.control.NumericEditField
        Label_2                        matlab.ui.control.Label
        TrasformarTypeDropDownLabel    matlab.ui.control.Label
        TrasformarTypeDropDown         matlab.ui.control.DropDown
        TrasformarMaterialDropDownLabel  matlab.ui.control.Label
        TrasformarMaterialDropDown     matlab.ui.control.DropDown
        ClearanceberweenHVYokecmEditFieldLabel  matlab.ui.control.Label
        ClearanceberweenHVYokecmEditField  matlab.ui.control.NumericEditField
        MaxAllowableTemperatureinCEditFieldLabel  matlab.ui.control.Label
        MaxAllowableTemperatureinCEditField  matlab.ui.control.NumericEditField
        OffLoadTapChangerEditFieldLabel  matlab.ui.control.Label
        OffLoadTapChangerEditField     matlab.ui.control.NumericEditField
        CoreTypeDropDownLabel          matlab.ui.control.Label
        CoreTypeDropDown               matlab.ui.control.DropDown
        ConnectionTypeDropDownLabel    matlab.ui.control.Label
        ConnectionTypeDropDown         matlab.ui.control.DropDown
        ExpectedPFEditFieldLabel       matlab.ui.control.Label
        ExpectedPFEditField            matlab.ui.control.NumericEditField
        ExpectedLoadingEditFieldLabel  matlab.ui.control.Label
        ExpectedLoadingEditField       matlab.ui.control.NumericEditField
        CoreConstructionDropDownLabel  matlab.ui.control.Label
        CoreConstructionDropDown       matlab.ui.control.DropDown
        FrequencyHzEditFieldLabel      matlab.ui.control.Label
        FrequencyHzEditField           matlab.ui.control.NumericEditField
        Label_3                        matlab.ui.control.Label
        Image                          matlab.ui.control.Image
        Label_4                        matlab.ui.control.Label
        StepUpTransformerDesignSimulatorLabel  matlab.ui.control.Label
        InputParametersLabel           matlab.ui.control.Label
        RightPanel                     matlab.ui.container.Panel
        BmTEditFieldLabel              matlab.ui.control.Label
        BmTEditField                   matlab.ui.control.NumericEditField
        SimulateDesignButton           matlab.ui.control.Button
        KEditFieldLabel                matlab.ui.control.Label
        KEditField                     matlab.ui.control.NumericEditField
        ConstantAssumedBasedOnYourtransformerTypeMaterialLabel  matlab.ui.control.Label
        CoreDimensionsResultsLabel     matlab.ui.control.Label
        WindingParametersLabel         matlab.ui.control.Label
        ReqEditFieldLabel              matlab.ui.control.Label
        ReqEditField                   matlab.ui.control.NumericEditField
        CircuitsElectricAndMagneticParametersLabel  matlab.ui.control.Label
        ZeqEditFieldLabel              matlab.ui.control.Label
        ZeqEditField                   matlab.ui.control.NumericEditField
        Zeq_puEditFieldLabel           matlab.ui.control.Label
        Zeq_puEditField                matlab.ui.control.NumericEditField
        Formar1010ConductorLabel       matlab.ui.control.Label
        LVWindingParametersLabel       matlab.ui.control.Label
        HVWindingParametersLabel       matlab.ui.control.Label
        HyokeAtmEditFieldLabel         matlab.ui.control.Label
        HyokeAtmEditField              matlab.ui.control.NumericEditField
        HlambAtmEditFieldLabel         matlab.ui.control.Label
        HlambAtmEditField              matlab.ui.control.NumericEditField
        VoltperturnEditFieldLabel      matlab.ui.control.Label
        VoltperturnEditField           matlab.ui.control.NumericEditField
        KwEditFieldLabel               matlab.ui.control.Label
        KwEditField                    matlab.ui.control.NumericEditField
        Aim2EditFieldLabel             matlab.ui.control.Label
        Aim2EditField                  matlab.ui.control.NumericEditField
        dcmEditFieldLabel              matlab.ui.control.Label
        dcmEditField                   matlab.ui.control.NumericEditField
        acmEditFieldLabel              matlab.ui.control.Label
        acmEditField                   matlab.ui.control.NumericEditField
        hyokecmEditFieldLabel          matlab.ui.control.Label
        hyokecmEditField               matlab.ui.control.NumericEditField
        WwcmEditFieldLabel             matlab.ui.control.Label
        WwcmEditField                  matlab.ui.control.NumericEditField
        HwcmEditFieldLabel             matlab.ui.control.Label
        HwcmEditField                  matlab.ui.control.NumericEditField
        DcmEditFieldLabel              matlab.ui.control.Label
        DcmEditField                   matlab.ui.control.NumericEditField
        CoreheighcmtEditFieldLabel     matlab.ui.control.Label
        CoreheighcmtEditField          matlab.ui.control.NumericEditField
        CoreWidthcmEditFieldLabel      matlab.ui.control.Label
        CoreWidthcmEditField           matlab.ui.control.NumericEditField
        Awm2EditFieldLabel             matlab.ui.control.Label
        Awm2EditField                  matlab.ui.control.NumericEditField
        CoolingMethodEditFieldLabel    matlab.ui.control.Label
        CoolingMethodEditField         matlab.ui.control.EditField
        HVturnsEditFieldLabel          matlab.ui.control.Label
        HVturnsEditField               matlab.ui.control.NumericEditField
        LVturnsEditFieldLabel          matlab.ui.control.Label
        LVturnsEditField               matlab.ui.control.NumericEditField
        HVcurrentAEditFieldLabel       matlab.ui.control.Label
        HVcurrentAEditField            matlab.ui.control.NumericEditField
        LVcurrentAEditFieldLabel       matlab.ui.control.Label
        LVcurrentAEditField            matlab.ui.control.NumericEditField
        HVconductorcsamm2EditFieldLabel  matlab.ui.control.Label
        HVconductorcsamm2EditField     matlab.ui.control.NumericEditField
        LVconductorcsamm2EditFieldLabel  matlab.ui.control.Label
        LVconductorcsamm2EditField     matlab.ui.control.NumericEditField
        LVinnerdiametercmEditFieldLabel  matlab.ui.control.Label
        LVinnerdiametercmEditField     matlab.ui.control.NumericEditField
        XeqEditFieldLabel              matlab.ui.control.Label
        XeqEditField                   matlab.ui.control.NumericEditField
        LVHeightcmEditFieldLabel       matlab.ui.control.Label
        LVHeightcmEditField            matlab.ui.control.NumericEditField
        FormarHeightcmEditFieldLabel   matlab.ui.control.Label
        FormarHeightcmEditField        matlab.ui.control.NumericEditField
        FormarwidthcmEditFieldLabel    matlab.ui.control.Label
        FormarwidthcmEditField         matlab.ui.control.NumericEditField
        HVHeightcmEditFieldLabel       matlab.ui.control.Label
        HVHeightcmEditField            matlab.ui.control.NumericEditField
        HVWidthcmEditFieldLabel        matlab.ui.control.Label
        HVWidthcmEditField             matlab.ui.control.NumericEditField
        ZbaseEditFieldLabel            matlab.ui.control.Label
        ZbaseEditField                 matlab.ui.control.NumericEditField
        HVconductortypeEditField_2Label  matlab.ui.control.Label
        HVconductortypeEditField_2     matlab.ui.control.EditField
        LVconductortypeEditFieldLabel  matlab.ui.control.Label
        LVconductortypeEditField       matlab.ui.control.EditField
        IcoreAEditFieldLabel           matlab.ui.control.Label
        IcoreAEditField                matlab.ui.control.NumericEditField
        InoloadAEditFieldLabel         matlab.ui.control.Label
        InoloadAEditField              matlab.ui.control.NumericEditField
        RcoreohmEditFieldLabel         matlab.ui.control.Label
        RcoreohmEditField              matlab.ui.control.NumericEditField
        ImagtAEditFieldLabel           matlab.ui.control.Label
        ImagtAEditField                matlab.ui.control.NumericEditField
        CopperlosswattEditFieldLabel   matlab.ui.control.Label
        CopperlosswattEditField        matlab.ui.control.NumericEditField
        ExpectedefficiencyEditFieldLabel  matlab.ui.control.Label
        ExpectedefficiencyEditField    matlab.ui.control.NumericEditField
        TankMeasurementsaDimensionsLabel  matlab.ui.control.Label
        TanklengthXmEditFieldLabel     matlab.ui.control.Label
        TanklengthXmEditField          matlab.ui.control.NumericEditField
        TankwidthYmEditFieldLabel      matlab.ui.control.Label
        TankwidthYmEditField           matlab.ui.control.NumericEditField
        TankheightZmEditFieldLabel     matlab.ui.control.Label
        TankheightZmEditField          matlab.ui.control.NumericEditField
        Tankdissipationareainm2EditFieldLabel  matlab.ui.control.Label
        Tankdissipationareainm2EditField  matlab.ui.control.NumericEditField
        MaximumtemperatureriseCLabel   matlab.ui.control.Label
        MaximumtemperatureriseCEditField  matlab.ui.control.NumericEditField
        CoolingtubenecessityEditFieldLabel  matlab.ui.control.Label
        CoolingtubenecessityEditField  matlab.ui.control.EditField
        DesignFeasibilityEditFieldLabel  matlab.ui.control.Label
        DesignFeasibilityEditField     matlab.ui.control.EditField
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
%     properties (Access = private)
%         s = app.TransformerRatingSinKVAEditField.Value*1000;
%         hv_volt=app.SecondaryVoltageHVinKVEditField.Value*1000;
%         lv_volt=app.PrimaryVoltageLVinVEditField.Value;
%         f=app.frquDropDown.Value;
%         cd=app.CurrentDensityinAmm2EditField.Value;
%         trans_type=app.TrasformarTypeDropDown.Value
%         trans_mat=app.TrasformarMaterialDropDown.Value;
%         taps=app.OffLoadTapChangerEditField.Value/100;
%         c=app.ClearanceberweenHVYokecmEditField.Value;
%         temp_max=app.MaxAllowableTemperatureinCEditField.Value;
%         core_type=app.CoreTypeDropDown.Value;
%         connection_type=app.ConnectionTypeDropDown.Value;
%         PF=app.ExpectedPFEditField.Value;
%         load_pre=app.ExpectedLoadingEditField.Value/100;
%     end
    
    methods (Access = private)
        
        function [Et,Ai,Kw,Aw,d,a,h_yoke,Ww,Hw,D,width,hight]=core_dim(app,s,hv_volt,f,bm,k,cd,core_construction)
           
            % First step using only s we get Et
            Et=k*sqrt(s/1000);
            %k is const  and S in KVA unit
            %2nd step Get A_iron using Et
            Ai = Et/(4.44*f*bm);
            %3nd S in VA in get Kw using HV side volt, AW using equ below
            
            %from S*1000=3.33*Kw*cd*Bm*A_iron*Aw
            if s <= 50000
                Kw = 8/((50+hv_volt/1000));
            elseif s <= 200000
                Kw = 10/(30+(hv_volt/1000));
            elseif s <= 1000000
                Kw = 12/(30+(hv_volt/1000));
            end
            Aw=(s)/(3.33*Kw*(cd*10^6)*bm*Ai*f);
            %cruciform core
            if contains(core_construction,'2 Stepped')
                d = sqrt(Ai*10^4/0.56);
            elseif contains(core_construction,'3 Stepped')
                d = sqrt(Ai*10^4/0.6);
            end    
            %4nd assume get d = sqrt(A_iron/.56) cruciform core
            a=.9*d;
            h_yoke=1.2*0.9*d;
            %5nd a=.9d,h_yoke=1.2a=1.2*09d in cm
            Ww=sqrt((Aw/3))*10^2;%cm`
            Hw=3*Ww;
            D=Ww+d;
            width=2*D+a;
            hight=Hw+2*h_yoke;
            
        end
        
        function [Bm, H_limb , H_yoke ,K] = assume_Bm(app,trans_type,trans_mat,hv_volt)
             
            if contains(trans_type,'Power') && (contains(trans_mat,'Hot Rolled'))
                Bm=1.35;
                H_limb=900;
                H_yoke=400;
                K= .65;
            elseif  contains(trans_type,'Distribution') && (contains(trans_mat,'Hot Rolled'))
                Bm=1.25 ;
                H_limb=550;
                H_yoke=285;
                K=.45;
            elseif  contains(trans_mat,'Cold Rolled')
                if hv_volt < (215*1000)
                    Bm=1.55;
                    H_limb=1125;
                    H_yoke=545;
                    K = 0.62;
                else
                    Bm=1.7;
                    H_limb=1500;
                    H_yoke=920;
                    K = 0.62;
                end
            end
        end
        
         % Description
        function cooling_ty = cooling_type(app,s)
            if   isnumeric(s)
                if s < 20000
                    cooling_ty='(AN) Cooling';
                elseif  s <= 50000
                    cooling_ty='(AB) Cooling';
                elseif  s<= 10000000
                    cooling_ty='(ON) Cooling';
                elseif  s<= 15000000
                    cooling_ty='(ONAF) Cooling';
                elseif  s <= 25000000
                    cooling_ty='(OFAN) Cooling';
                elseif s <= 60000000
                    cooling_ty='(OFAF) Cooling';
                elseif s > 60000000
                    cooling_ty='(OFWF) Cooling';
                end
            end
              %             if (1.1>=min_air)&&(cd >max_air)
%                 cooling_ty='Air Natural';
%             elseif (min_oil>=cd)&&(max_oil >3.2)
%                 cooling_ty='Oil Natural';
%                 
%             elseif cd>max_oil
%                 cooling_ty ='forced';
%             end
        end
        
        
       
        
        function [hv_volt_ph,lv_volt_ph] = phase_volt(app,connection_type,hv_volt,lv_volt)
            switch connection_type
                case 'Deltaÿ,StarY'
                     hv_volt_ph = hv_volt;
                     lv_volt_ph = lv_volt/sqrt(3);
                case 'Deltaÿ,Deltaÿ'
                    hv_volt_ph = hv_volt;
                    lv_volt_ph = lv_volt;
                case 'StarY,Deltaÿ'
                    hv_volt_ph = hv_volt/sqrt(3);
                    lv_volt_ph = lv_volt;
                    
                case  'StarY,StarY'  
                    hv_volt_ph = hv_volt/sqrt(3);
                    lv_volt_ph = lv_volt/sqrt(3);
            end
        end
        
        
        function [Req_ref_hv,Xeq_ref_hv,Z_eq,Z_base,Z_eq_pu,turn_ratio] = Z_eq_calc(app,d,d_outer_lv,T_l,T_h,c_csa_lv,c_csa_hv,c_w_lv,d_t_hv,hv_volt,I_hv,d_in_hv,d_out_hv,f)
                 L_mean_lv=(pi*(((d+2)+d_outer_lv)/2))*10^-2;%m
                 R_lv=.021*((T_l*L_mean_lv)/c_csa_lv) ;
                 
                 L_mean_hv=(pi*((d_in_hv+2+d_out_hv)/2))*10^-2;%m
                 R_hv=.021*((T_h*L_mean_hv)/(c_csa_hv)) ;%ohm
                 
                 
                 %Req cac
                 turn_ratio=T_h/T_l;
                 Req_ref_hv=R_hv+(R_lv*(turn_ratio^2));
                 %Xeq calc
                 L_c=(L_mean_lv+L_mean_hv)/2;
                 L_mean_hv_lv = pi*(((d+2+d_out_hv)*10^-2)/2);
                 
                 b1= (c_w_lv+4)*10^-3;
                 b2= d_t_hv*10*10^-3;
                 clearance= .06;
                 
                 Xeq_ref_hv=2*pi*f*4*pi*(10^-7)*(T_h^2)*(L_mean_hv_lv/L_c)*(clearance+((b1+b2)/3));
                 %Zeq calc
                 Z_eq=sqrt(((Req_ref_hv^2)+(Xeq_ref_hv^2)));
                 Z_base=(hv_volt/I_hv);
                 Z_eq_pu=Z_eq/Z_base;
            
        end
        
        
        
        
        function [T_h,T_l,T_tap,I_hv,I_lv,c_csa_hv,c_csa_lv,conductor_type_hv,conductor_type_lv,d_outer_lv,d_outer_f,Lc_hv_10,Lc_tot_hv,d_out_hv,Lc_lv,c_w_lv,d_t_hv,d_in_hv] = winding_parameters(app,Et,s,hv_volt,hv_volt_ph,lv_volt,lv_volt_ph,taps,cd,d)
            T_h = round((hv_volt_ph)/(Et)) ;
            T_tap =round( T_h * (taps));
            T_l= round((lv_volt_ph)/(Et)) ;
            
            I_hv=(s)/(3*hv_volt_ph);
            
            I_lv=(s)/(3*lv_volt_ph);
            
            c_csa_hv=I_hv/cd; %in mm2
            c_csa_lv=I_lv/cd; % in mm2
            
            conductor_type_hv ='Circular'; % case of 10*10 formar H-V side
            % c_csa_s=pi/4 * d^2
            d_t_hv=sqrt((4*c_csa_hv)/pi)+2;%mm
            Lc_hv_10=10*(d_t_hv*10^-1)+2;%cm
            d_outer_f = Lc_hv_10-1;%Formar Width
            % end
            
            
            %if c_csa_p > 10 L-V side
            conductor_type_lv ='Rectangular ';%4*3 rectanlar
            c_w_lv = ceil(sqrt(c_csa_lv));%mm%condactor width
            c_h_lv = floor(sqrt(c_csa_lv));%mm%condactor Hight
            d_outer_lv = d+2+(((c_w_lv+4)*10^-1)*2);%cm
            Lc_lv = (c_h_lv+5)*(10^-1)*T_l;%cm
            
            %
            d_in_hv = d_outer_lv+10;%cm
            d_out_hv = d_in_hv+(d_outer_f*2);%cm
            %check d_outer_p = d+2+((c_w_lv+4)*10^-1)<HW if not redesign
            % get height of H-V side wich include sum of fromar height
            t_h =T_h+T_tap;
            Lc_tot_hv =0;%sum of formar height == height of H-V side 
            
            while t_h >= 10
                if t_h>= 100
                    Lc_tot_hv=Lc_tot_hv+Lc_hv_10;
                    t_h=t_h-100;
                else
                    if rem(t_h , 10 ) == 0
                        Lc_tot_hv=Lc_tot_hv+((t_h/10)*(d_t_hv*10^-1)+2);
                        t_h=0;
                    else
                        x= (t_h-rem(t_h , 10 ))/10;
                        Lc_tot_hv=Lc_tot_hv+((x*(d_t_hv*10^-1)+2)) ;
                        t_h= rem(t_h , 10 );
                    end
                end
                
            end
            % winding_parameters
           
            
            
        end
        
       
        
        
        
        function [PI_tot,R_core,I_core,I_mag,I_nl,eff,P_cu_fl] = Loses_eff(app,s,Bm,Ai,Hw,width,H_limb,hv_volt_ph,H_yoke,T_h,PF,load_pre,Req_ref_hv)
            %core loss Current calc
            %limb calc
            Bm_limb=Bm;
            Pi_limb=2.6;%watt/kg
            vol_limb=3*Ai*(Hw*10^-2);%m3
            mass_limb=7800*vol_limb;%kg
            %yoke calc
            Bm_yoke = Bm*.8;
            Pi_yoke = 1.3;
            vol_yoke=2*(width*10^-2)*(Ai*1.2);
            mass_yoke=7800*vol_yoke;
            %for 3-ph
            PI_tot=(Pi_limb*mass_limb)+(mass_yoke*Pi_yoke);
            
            R_core=((3*(hv_volt_ph^2))/PI_tot);
            
            I_core=hv_volt_ph/R_core;
            
            %magnetizing Current calc
            AT=((3*Hw*(10^-2)*H_limb)+(2*H_yoke*(width*(10^-2)))+(Bm_limb/(4*pi*10^-7))*6*.5*10^-3);
            %Im calc
            I_mag=(AT/(3*sqrt(2)*T_h));
            
            %calc No-Load Current
            I_nl=sqrt(((I_mag^2)+(I_core^2)));
            %regulation and Eff calc
            P_cu_fl=3*((s/(3*hv_volt_ph))^2)*Req_ref_hv;
            
            eff =(load_pre*s*PF)/((load_pre*s*PF)+(PI_tot)+(load_pre^2*P_cu_fl));
            
            %regu=x*(((Req_ref_hv/Z_base)*PF)+((xeq_ref_hv/Z_eq)*sin(acos(PF))));
        end
        
        
        
        function [X,Y,Z,tank_diss_area,temp_rise_max,tube] = tank_dim(app,d_out_hv,width,hight,a,temp_max,PI_tot,P_cu_fl,hv_volt)
            %tank design
           
            if hv_volt <= 66000
            Dt = 12;
            else
            Dt = 32;
            end
            
            X=(d_out_hv+(2*Dt))*10^-2;%m
            Y=(width+d_out_hv+(2*Dt)-a)*10^-2;
            Z=(hight*10^-2)+.5;
            tank_diss_area=(2*X*Z)+(2*Y*Z)+(X*Y)/2;
            
            %maximum temperatureMax temp
            temp_rise_max=(PI_tot+P_cu_fl)/(12.5*tank_diss_area);
            if  temp_max < temp_rise_max%temp_max enterd by the user
                tube='cooling tubes is needed';
            else
                tube=' No need for cooling tubes';
            end
           
            
        end
        
        
        
        
        
        
        
        function results = suggest_redesign(app,Hw,Lc_lv,Lc_tot_hv,c,D,d_out_hv,Z_eq_pu,I_nl,I_hv)
            if (Lc_lv < Hw) && ((Lc_tot_hv + 2*c) < Hw ) && (D >= (d_out_hv +2*c)) && (Z_eq_pu <= 0.1) && (I_nl/I_hv <= 0.05)
                results = 'Not Feasible Design';
            else
                results ='Feasible Design';
                
            end
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SimulateDesignButton
        function SimulateDesignButtonPushed(app, event)
            s = app.TransformerRatingSinKVAEditField.Value*1000;
            hv_volt=app.SecondaryVoltageHVinKVEditField.Value*1000;
            lv_volt=app.PrimaryVoltageLVinVEditField.Value;
            f=app.FrequencyHzEditField.Value;
            cd=app.CurrentDensityinAmm2EditField.Value;
            trans_type=app.TrasformarTypeDropDown.Value;
            trans_mat=app.TrasformarMaterialDropDown.Value;
            taps=app.OffLoadTapChangerEditField.Value/100;
            c=app.ClearanceberweenHVYokecmEditField.Value;
            temp_max=app.MaxAllowableTemperatureinCEditField.Value;
            core_type=app.CoreTypeDropDown.Value;
            connection_type=app.ConnectionTypeDropDown.Value;
            PF=app.ExpectedPFEditField.Value;
            load_pre=app.ExpectedLoadingEditField.Value/100;
            core_construction=app.CoreConstructionDropDown.Value;
            %-------------------------------------------------------------
            %data prepration
            [hv_volt_ph , lv_volt_ph] = phase_volt(app,connection_type,hv_volt,lv_volt);
            %----------------------------------consttant------------------
            [bm, H_limb , H_yoke ,K] = assume_Bm(app,trans_type,trans_mat,hv_volt);
%             bm= r(1,1);
%             H_limb = r(1,2) ;
%             H_yoke =r(1,3);
%             k=r(1,4);
            app.HyokeAtmEditField.Value = H_yoke;
            app.HlambAtmEditField.Value = H_limb;
            app.KEditField.Value = K;
            app.BmTEditField.Value = bm;
            %-----------------------------------------------------------
            cooling_ty = cooling_type(app,s);
            app.CoolingMethodEditField.Value = cooling_ty; 
            %------------------------------------------------------------
            
            [Et,Ai,Kw,Aw,d,a,h_yoke,Ww,Hw,D,width,hight]=core_dim(app,s,hv_volt,f,bm,K,cd,core_construction);
            
            app.VoltperturnEditField.Value = Et;
            
            app.Aim2EditField.Value = Ai;
            
            app.KwEditField.Value = Kw;
            
            app.Awm2EditField.Value = Aw;
            
            app.dcmEditField.Value = d;
            
            app.acmEditField.Value = a;
            
            app.hyokecmEditField.Value = h_yoke;
            
            app.WwcmEditField.Value = Ww;
            
            app.HwcmEditField.Value =Hw;
            
            app.DcmEditField.Value = D;
            
            app.CoreWidthcmEditField.Value = width;
            
            app.CoreheighcmtEditField.Value = hight;
%--------------------------------------------------------------------------            
            [T_h,T_l,T_tap,I_hv,I_lv,c_csa_hv,c_csa_lv,conductor_type_hv,conductor_type_lv,d_outer_lv,d_outer_f,Lc_hv_10,Lc_tot_hv,d_out_hv,Lc_lv,c_w_lv,d_t_hv,d_in_hv]= winding_parameters(app,Et,s,hv_volt,hv_volt_ph,lv_volt,lv_volt_ph,taps,cd,d);

            app.HVturnsEditField.Value = T_h;
            app.LVturnsEditField.Value = T_l;
            app.HVcurrentAEditField.Value = I_hv;
            app.LVcurrentAEditField.Value = I_lv;
            app.HVconductorcsamm2EditField.Value = c_csa_hv;
            app.LVconductorcsamm2EditField.Value = c_csa_lv;
            %
            app.HVconductortypeEditField_2.Value = conductor_type_hv;
            app.LVconductortypeEditField.Value = conductor_type_lv;
            app.LVinnerdiametercmEditField.Value = d_outer_lv;
            app.HVHeightcmEditField.Value = Lc_tot_hv;
            app.HVWidthcmEditField.Value = d_out_hv;
            app.LVHeightcmEditField.Value = Lc_lv;
            %
            app.FormarHeightcmEditField.Value = Lc_hv_10;
            app.FormarwidthcmEditField.Value = d_outer_f;
%-------------------------------------------------------------------------
            [Req_ref_hv,Xeq_ref_hv,Z_eq,Z_base,Z_eq_pu,turn_ratio] = Z_eq_calc(app,d,d_outer_lv,T_l,T_h,c_csa_lv,c_csa_hv,c_w_lv,d_t_hv,hv_volt,I_hv,d_in_hv,d_out_hv,f);
            
            app.ReqEditField.Value = Req_ref_hv;
            app.XeqEditField.Value = Xeq_ref_hv;
            %
            app.ZeqEditField.Value = Z_eq;
            app.ZbaseEditField.Value = Z_base;
            app.Zeq_puEditField.Value = Z_eq_pu;
%---------------------------------------------------------------------------------------------------------------------------            
            [PI_tot,R_core,I_core,I_mag,I_nl,eff,P_cu_fl] = Loses_eff(app,s,bm,Ai,Hw,width,H_limb,hv_volt_ph,H_yoke,T_h,PF,load_pre,Req_ref_hv);
            
            app.IcoreAEditField.Value = I_core;
            app.InoloadAEditField.Value = I_nl;
            app.RcoreohmEditField.Value = R_core;
            app.ImagtAEditField.Value = I_mag;
            app.CopperlosswattEditField.Value = P_cu_fl;
            app.ExpectedefficiencyEditField.Value = eff*100;
%-------------------------------------------------------------------------------------------------------------------------------            
            [X,Y,Z,tank_diss_area,temp_rise_max,tube] = tank_dim(app,d_out_hv,width,hight,a,temp_max,PI_tot,P_cu_fl,hv_volt);
            
            app.TanklengthXmEditField.Value = X;
            app.TankwidthYmEditField.Value = Y;
            app.TankheightZmEditField.Value = Z;
            app.Tankdissipationareainm2EditField.Value = tank_diss_area;
            app.MaximumtemperatureriseCEditField.Value = temp_rise_max;
            app.CoolingtubenecessityEditField.Value = tube;
%--------------------------------------------------------------------------------------------------------------------------------            
            
            
            results = suggest_redesign(app,Hw,Lc_lv,Lc_tot_hv,c,D,d_out_hv,Z_eq_pu,I_nl,I_hv);
            app.DesignFeasibilityEditField.Value =results;
                
            
            
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {964, 964};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {502, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1103 964];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {502, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create TransformerRatingSinKVAEditFieldLabel
            app.TransformerRatingSinKVAEditFieldLabel = uilabel(app.LeftPanel);
            app.TransformerRatingSinKVAEditFieldLabel.HorizontalAlignment = 'right';
            app.TransformerRatingSinKVAEditFieldLabel.Position = [16 492 182 22];
            app.TransformerRatingSinKVAEditFieldLabel.Text = 'Transformer Rating  ''S'' in (KVA) ';

            % Create TransformerRatingSinKVAEditField
            app.TransformerRatingSinKVAEditField = uieditfield(app.LeftPanel, 'numeric');
            app.TransformerRatingSinKVAEditField.Limits = [0 Inf];
            app.TransformerRatingSinKVAEditField.Position = [216 492 266 22];

            % Create PrimaryVoltageLVinVLabel
            app.PrimaryVoltageLVinVLabel = uilabel(app.LeftPanel);
            app.PrimaryVoltageLVinVLabel.HorizontalAlignment = 'right';
            app.PrimaryVoltageLVinVLabel.Position = [24 454 149 22];
            app.PrimaryVoltageLVinVLabel.Text = 'Primary Voltage ''L-V'' in (V)';

            % Create PrimaryVoltageLVinVEditField
            app.PrimaryVoltageLVinVEditField = uieditfield(app.LeftPanel, 'numeric');
            app.PrimaryVoltageLVinVEditField.Position = [216 454 263 22];

            % Create SecondaryVoltageHVinKVLabel
            app.SecondaryVoltageHVinKVLabel = uilabel(app.LeftPanel);
            app.SecondaryVoltageHVinKVLabel.HorizontalAlignment = 'right';
            app.SecondaryVoltageHVinKVLabel.Position = [16 412 170 22];
            app.SecondaryVoltageHVinKVLabel.Text = 'Secondary Voltage ''H-V)''in (KV)';

            % Create SecondaryVoltageHVinKVEditField
            app.SecondaryVoltageHVinKVEditField = uieditfield(app.LeftPanel, 'numeric');
            app.SecondaryVoltageHVinKVEditField.Position = [220 412 266 22];

            % Create CurrentDensityinAmm2EditFieldLabel
            app.CurrentDensityinAmm2EditFieldLabel = uilabel(app.LeftPanel);
            app.CurrentDensityinAmm2EditFieldLabel.HorizontalAlignment = 'right';
            app.CurrentDensityinAmm2EditFieldLabel.Position = [32 173 151 22];
            app.CurrentDensityinAmm2EditFieldLabel.Text = 'Current Density in (A/mm2)';

            % Create CurrentDensityinAmm2EditField
            app.CurrentDensityinAmm2EditField = uieditfield(app.LeftPanel, 'numeric');
            app.CurrentDensityinAmm2EditField.Position = [217 173 266 22];

            % Create Label_2
            app.Label_2 = uilabel(app.LeftPanel);
            app.Label_2.Position = [458 725 25 22];
            app.Label_2.Text = '';

            % Create TrasformarTypeDropDownLabel
            app.TrasformarTypeDropDownLabel = uilabel(app.LeftPanel);
            app.TrasformarTypeDropDownLabel.HorizontalAlignment = 'right';
            app.TrasformarTypeDropDownLabel.Position = [60 331 94 22];
            app.TrasformarTypeDropDownLabel.Text = 'Trasformar Type';

            % Create TrasformarTypeDropDown
            app.TrasformarTypeDropDown = uidropdown(app.LeftPanel);
            app.TrasformarTypeDropDown.Items = {'Power', 'Distribution'};
            app.TrasformarTypeDropDown.Position = [220 331 266 22];
            app.TrasformarTypeDropDown.Value = 'Power';

            % Create TrasformarMaterialDropDownLabel
            app.TrasformarMaterialDropDownLabel = uilabel(app.LeftPanel);
            app.TrasformarMaterialDropDownLabel.HorizontalAlignment = 'right';
            app.TrasformarMaterialDropDownLabel.Position = [57 364 105 22];
            app.TrasformarMaterialDropDownLabel.Text = 'Trasformar Material';

            % Create TrasformarMaterialDropDown
            app.TrasformarMaterialDropDown = uidropdown(app.LeftPanel);
            app.TrasformarMaterialDropDown.Items = {'Hot Rolled', 'Cold Rolled'};
            app.TrasformarMaterialDropDown.Position = [219 364 266 22];
            app.TrasformarMaterialDropDown.Value = 'Hot Rolled';

            % Create ClearanceberweenHVYokecmEditFieldLabel
            app.ClearanceberweenHVYokecmEditFieldLabel = uilabel(app.LeftPanel);
            app.ClearanceberweenHVYokecmEditFieldLabel.HorizontalAlignment = 'right';
            app.ClearanceberweenHVYokecmEditFieldLabel.Position = [-3 144 195 22];
            app.ClearanceberweenHVYokecmEditFieldLabel.Text = 'Clearance berween HV & Yoke(cm)';

            % Create ClearanceberweenHVYokecmEditField
            app.ClearanceberweenHVYokecmEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ClearanceberweenHVYokecmEditField.Position = [215 144 266 22];

            % Create MaxAllowableTemperatureinCEditFieldLabel
            app.MaxAllowableTemperatureinCEditFieldLabel = uilabel(app.LeftPanel);
            app.MaxAllowableTemperatureinCEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxAllowableTemperatureinCEditFieldLabel.Position = [21 16 182 22];
            app.MaxAllowableTemperatureinCEditFieldLabel.Text = 'Max Allowable Temperature in °C';

            % Create MaxAllowableTemperatureinCEditField
            app.MaxAllowableTemperatureinCEditField = uieditfield(app.LeftPanel, 'numeric');
            app.MaxAllowableTemperatureinCEditField.Position = [219 16 266 22];

            % Create OffLoadTapChangerEditFieldLabel
            app.OffLoadTapChangerEditFieldLabel = uilabel(app.LeftPanel);
            app.OffLoadTapChangerEditFieldLabel.HorizontalAlignment = 'right';
            app.OffLoadTapChangerEditFieldLabel.Position = [38 114 138 22];
            app.OffLoadTapChangerEditFieldLabel.Text = 'Off-Load Tap Changer %';

            % Create OffLoadTapChangerEditField
            app.OffLoadTapChangerEditField = uieditfield(app.LeftPanel, 'numeric');
            app.OffLoadTapChangerEditField.Position = [219 114 267 22];

            % Create CoreTypeDropDownLabel
            app.CoreTypeDropDownLabel = uilabel(app.LeftPanel);
            app.CoreTypeDropDownLabel.HorizontalAlignment = 'right';
            app.CoreTypeDropDownLabel.Position = [86 269 61 22];
            app.CoreTypeDropDownLabel.Text = 'Core Type';

            % Create CoreTypeDropDown
            app.CoreTypeDropDown = uidropdown(app.LeftPanel);
            app.CoreTypeDropDown.Items = {'Core', 'Shell '};
            app.CoreTypeDropDown.Position = [224 269 266 22];
            app.CoreTypeDropDown.Value = 'Core';

            % Create ConnectionTypeDropDownLabel
            app.ConnectionTypeDropDownLabel = uilabel(app.LeftPanel);
            app.ConnectionTypeDropDownLabel.HorizontalAlignment = 'right';
            app.ConnectionTypeDropDownLabel.Position = [68 299 93 22];
            app.ConnectionTypeDropDownLabel.Text = 'Connection Type';

            % Create ConnectionTypeDropDown
            app.ConnectionTypeDropDown = uidropdown(app.LeftPanel);
            app.ConnectionTypeDropDown.Items = {'Deltaÿ,StarY', 'Deltaÿ,Deltaÿ', 'StarY,Deltaÿ', 'StarY,StarY'};
            app.ConnectionTypeDropDown.Position = [224 299 266 22];
            app.ConnectionTypeDropDown.Value = 'Deltaÿ,StarY';

            % Create ExpectedPFEditFieldLabel
            app.ExpectedPFEditFieldLabel = uilabel(app.LeftPanel);
            app.ExpectedPFEditFieldLabel.HorizontalAlignment = 'right';
            app.ExpectedPFEditFieldLabel.Position = [70 82 74 22];
            app.ExpectedPFEditFieldLabel.Text = 'Expected PF';

            % Create ExpectedPFEditField
            app.ExpectedPFEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ExpectedPFEditField.Position = [222 82 266 22];

            % Create ExpectedLoadingEditFieldLabel
            app.ExpectedLoadingEditFieldLabel = uilabel(app.LeftPanel);
            app.ExpectedLoadingEditFieldLabel.HorizontalAlignment = 'right';
            app.ExpectedLoadingEditFieldLabel.Position = [50 47 115 22];
            app.ExpectedLoadingEditFieldLabel.Text = 'Expected Loading %';

            % Create ExpectedLoadingEditField
            app.ExpectedLoadingEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ExpectedLoadingEditField.Position = [219 47 266 22];

            % Create CoreConstructionDropDownLabel
            app.CoreConstructionDropDownLabel = uilabel(app.LeftPanel);
            app.CoreConstructionDropDownLabel.HorizontalAlignment = 'right';
            app.CoreConstructionDropDownLabel.Position = [60 238 102 22];
            app.CoreConstructionDropDownLabel.Text = 'Core Construction';

            % Create CoreConstructionDropDown
            app.CoreConstructionDropDown = uidropdown(app.LeftPanel);
            app.CoreConstructionDropDown.Items = {'2 Stepped', '3 Stepped'};
            app.CoreConstructionDropDown.Position = [224 238 266 22];
            app.CoreConstructionDropDown.Value = '2 Stepped';

            % Create FrequencyHzEditFieldLabel
            app.FrequencyHzEditFieldLabel = uilabel(app.LeftPanel);
            app.FrequencyHzEditFieldLabel.HorizontalAlignment = 'right';
            app.FrequencyHzEditFieldLabel.Position = [65 208 85 22];
            app.FrequencyHzEditFieldLabel.Text = 'Frequency(Hz)';

            % Create FrequencyHzEditField
            app.FrequencyHzEditField = uieditfield(app.LeftPanel, 'numeric');
            app.FrequencyHzEditField.Position = [218 208 267 22];

            % Create Label_3
            app.Label_3 = uilabel(app.LeftPanel);
            app.Label_3.BackgroundColor = [1 1 1];
            app.Label_3.FontSize = 11;
            app.Label_3.Position = [6 839 363 120];
            app.Label_3.Text = {'AIN SHAMS UNIVERSITY'; 'FACULTY OF ENGINEERING'; 'Department of Electrical Power and Machines, 3nd Year, Electrical'; 'Engineering '; ''; 'Research Project                            Course Code: EPM 322'; ''; 'Submitted by:ÿ'; 'Mostafa Mohamed Abdl-Aziz Badran'; 'ID:1601430                Bn:274ÿ'; ''};

            % Create Image
            app.Image = uiimage(app.LeftPanel);
            app.Image.BackgroundColor = [1 1 1];
            app.Image.Position = [367 840 129 119];
            app.Image.ImageSource = 'school_logo.jpg';

            % Create Label_4
            app.Label_4 = uilabel(app.LeftPanel);
            app.Label_4.FontSize = 15;
            app.Label_4.Position = [9 578 510 204];
            app.Label_4.Text = {'This project is part of course (Machines (2) EPM 322) and it purpose is to'; ' compute ÿElectric Transformer design based on the knowledge gain which '; 'had been taught by Prof. Adel Yousef Hannalla, Dr. Nabil Mohamed '; 'Hamed Mohamed,Dr. Ahmed Mohy ÿEldeen Ibraheem.ÿ This computational'; 'system is aiming to minimizing results error which can happenWhile doing'; 'design step  manually and notify with design feasibility and more features '; 'ÿitÿs core functions have been built by MATLAB language and for Graphical'; 'user interface(GUI) I used MATLAB app designer, but this system produce'; ' it result under some assumption ÿ'; 'ÿ1- it Step-Up Transformer ÿ'; 'ÿ2-tank design is plain walled tank shape'; '3-ÿTransformer electric circuit is circuit approximate and refer to H-V side ÿ'};

            % Create StepUpTransformerDesignSimulatorLabel
            app.StepUpTransformerDesignSimulatorLabel = uilabel(app.LeftPanel);
            app.StepUpTransformerDesignSimulatorLabel.HorizontalAlignment = 'center';
            app.StepUpTransformerDesignSimulatorLabel.FontSize = 20;
            app.StepUpTransformerDesignSimulatorLabel.FontWeight = 'bold';
            app.StepUpTransformerDesignSimulatorLabel.Position = [51 805 394 22];
            app.StepUpTransformerDesignSimulatorLabel.Text = ' Step-Up Transformer Design Simulator';

            % Create InputParametersLabel
            app.InputParametersLabel = uilabel(app.LeftPanel);
            app.InputParametersLabel.HorizontalAlignment = 'center';
            app.InputParametersLabel.FontWeight = 'bold';
            app.InputParametersLabel.Position = [30 546 464 22];
            app.InputParametersLabel.Text = '-------------------------------------------Input Parameters---------------------------------------------';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create BmTEditFieldLabel
            app.BmTEditFieldLabel = uilabel(app.RightPanel);
            app.BmTEditFieldLabel.HorizontalAlignment = 'right';
            app.BmTEditFieldLabel.Position = [14 874 39 22];
            app.BmTEditFieldLabel.Text = 'Bm(T)';

            % Create BmTEditField
            app.BmTEditField = uieditfield(app.RightPanel, 'numeric');
            app.BmTEditField.Position = [149 874 115 22];

            % Create SimulateDesignButton
            app.SimulateDesignButton = uibutton(app.RightPanel, 'push');
            app.SimulateDesignButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateDesignButtonPushed, true);
            app.SimulateDesignButton.FontSize = 15;
            app.SimulateDesignButton.FontWeight = 'bold';
            app.SimulateDesignButton.Position = [149 5 306 26];
            app.SimulateDesignButton.Text = 'Simulate Design';

            % Create KEditFieldLabel
            app.KEditFieldLabel = uilabel(app.RightPanel);
            app.KEditFieldLabel.HorizontalAlignment = 'right';
            app.KEditFieldLabel.Position = [18 847 11 22];
            app.KEditFieldLabel.Text = 'K';

            % Create KEditField
            app.KEditField = uieditfield(app.RightPanel, 'numeric');
            app.KEditField.Position = [150 847 114 22];

            % Create ConstantAssumedBasedOnYourtransformerTypeMaterialLabel
            app.ConstantAssumedBasedOnYourtransformerTypeMaterialLabel = uilabel(app.RightPanel);
            app.ConstantAssumedBasedOnYourtransformerTypeMaterialLabel.FontWeight = 'bold';
            app.ConstantAssumedBasedOnYourtransformerTypeMaterialLabel.Position = [10 895 587 22];
            app.ConstantAssumedBasedOnYourtransformerTypeMaterialLabel.Text = '.............................. Constant Assumed Based On Your transformer Type  & Material ..............................';

            % Create CoreDimensionsResultsLabel
            app.CoreDimensionsResultsLabel = uilabel(app.RightPanel);
            app.CoreDimensionsResultsLabel.FontWeight = 'bold';
            app.CoreDimensionsResultsLabel.FontColor = [0 0.4471 0.7412];
            app.CoreDimensionsResultsLabel.Position = [26 819 555 22];
            app.CoreDimensionsResultsLabel.Text = '.......................................................... Core Dimensions Results ..........................................................';

            % Create WindingParametersLabel
            app.WindingParametersLabel = uilabel(app.RightPanel);
            app.WindingParametersLabel.FontName = 'Monospaced';
            app.WindingParametersLabel.FontWeight = 'bold';
            app.WindingParametersLabel.FontColor = [1 0 0];
            app.WindingParametersLabel.Position = [40 644 525 22];
            app.WindingParametersLabel.Text = '.......................................................... Winding Parameters ..........................................................';

            % Create ReqEditFieldLabel
            app.ReqEditFieldLabel = uilabel(app.RightPanel);
            app.ReqEditFieldLabel.HorizontalAlignment = 'right';
            app.ReqEditFieldLabel.Position = [17 289 36 22];
            app.ReqEditFieldLabel.Text = 'Req';

            % Create ReqEditField
            app.ReqEditField = uieditfield(app.RightPanel, 'numeric');
            app.ReqEditField.Position = [65 289 62 22];

            % Create CircuitsElectricAndMagneticParametersLabel
            app.CircuitsElectricAndMagneticParametersLabel = uilabel(app.RightPanel);
            app.CircuitsElectricAndMagneticParametersLabel.FontWeight = 'bold';
            app.CircuitsElectricAndMagneticParametersLabel.FontColor = [0.4667 0.6745 0.1882];
            app.CircuitsElectricAndMagneticParametersLabel.Position = [15 331 575 22];
            app.CircuitsElectricAndMagneticParametersLabel.Text = '.............................................. Circuit''s Electric And Magnetic Parameters ..............................................';

            % Create ZeqEditFieldLabel
            app.ZeqEditFieldLabel = uilabel(app.RightPanel);
            app.ZeqEditFieldLabel.HorizontalAlignment = 'right';
            app.ZeqEditFieldLabel.Position = [287 290 21 22];
            app.ZeqEditFieldLabel.Text = 'Zeq';

            % Create ZeqEditField
            app.ZeqEditField = uieditfield(app.RightPanel, 'numeric');
            app.ZeqEditField.Position = [336 290 65 22];

            % Create Zeq_puEditFieldLabel
            app.Zeq_puEditFieldLabel = uilabel(app.RightPanel);
            app.Zeq_puEditFieldLabel.HorizontalAlignment = 'right';
            app.Zeq_puEditFieldLabel.Position = [148 289 46 22];
            app.Zeq_puEditFieldLabel.Text = 'Zeq_pu';

            % Create Zeq_puEditField
            app.Zeq_puEditField = uieditfield(app.RightPanel, 'numeric');
            app.Zeq_puEditField.Position = [209 290 55 22];

            % Create Formar1010ConductorLabel
            app.Formar1010ConductorLabel = uilabel(app.RightPanel);
            app.Formar1010ConductorLabel.FontSize = 11;
            app.Formar1010ConductorLabel.FontWeight = 'bold';
            app.Formar1010ConductorLabel.FontColor = [1 0 0];
            app.Formar1010ConductorLabel.Position = [10 385 134 22];
            app.Formar1010ConductorLabel.Text = 'Formar 10*10 Conductor';

            % Create LVWindingParametersLabel
            app.LVWindingParametersLabel = uilabel(app.RightPanel);
            app.LVWindingParametersLabel.FontSize = 11;
            app.LVWindingParametersLabel.FontWeight = 'bold';
            app.LVWindingParametersLabel.FontColor = [1 0 0];
            app.LVWindingParametersLabel.Position = [10 480 128 22];
            app.LVWindingParametersLabel.Text = 'LV Winding Parameters';

            % Create HVWindingParametersLabel
            app.HVWindingParametersLabel = uilabel(app.RightPanel);
            app.HVWindingParametersLabel.FontSize = 11;
            app.HVWindingParametersLabel.FontWeight = 'bold';
            app.HVWindingParametersLabel.FontColor = [1 0 0];
            app.HVWindingParametersLabel.Position = [10 608 134 22];
            app.HVWindingParametersLabel.Text = 'H-V Winding Parameters';

            % Create HyokeAtmEditFieldLabel
            app.HyokeAtmEditFieldLabel = uilabel(app.RightPanel);
            app.HyokeAtmEditFieldLabel.HorizontalAlignment = 'right';
            app.HyokeAtmEditFieldLabel.Position = [282 847 71 22];
            app.HyokeAtmEditFieldLabel.Text = 'H yoke(At/m)';

            % Create HyokeAtmEditField
            app.HyokeAtmEditField = uieditfield(app.RightPanel, 'numeric');
            app.HyokeAtmEditField.Position = [465 847 116 22];

            % Create HlambAtmEditFieldLabel
            app.HlambAtmEditFieldLabel = uilabel(app.RightPanel);
            app.HlambAtmEditFieldLabel.HorizontalAlignment = 'right';
            app.HlambAtmEditFieldLabel.Position = [283 874 72 22];
            app.HlambAtmEditFieldLabel.Text = 'H lamb(At/m)';

            % Create HlambAtmEditField
            app.HlambAtmEditField = uieditfield(app.RightPanel, 'numeric');
            app.HlambAtmEditField.Position = [466 874 116 22];

            % Create VoltperturnEditFieldLabel
            app.VoltperturnEditFieldLabel = uilabel(app.RightPanel);
            app.VoltperturnEditFieldLabel.HorizontalAlignment = 'right';
            app.VoltperturnEditFieldLabel.Position = [393 797 70 22];
            app.VoltperturnEditFieldLabel.Text = 'Volt per turn';

            % Create VoltperturnEditField
            app.VoltperturnEditField = uieditfield(app.RightPanel, 'numeric');
            app.VoltperturnEditField.Position = [496 797 85 22];

            % Create KwEditFieldLabel
            app.KwEditFieldLabel = uilabel(app.RightPanel);
            app.KwEditFieldLabel.HorizontalAlignment = 'right';
            app.KwEditFieldLabel.Position = [209 710 25 22];
            app.KwEditFieldLabel.Text = 'Kw';

            % Create KwEditField
            app.KwEditField = uieditfield(app.RightPanel, 'numeric');
            app.KwEditField.Position = [302 710 84 22];

            % Create Aim2EditFieldLabel
            app.Aim2EditFieldLabel = uilabel(app.RightPanel);
            app.Aim2EditFieldLabel.HorizontalAlignment = 'right';
            app.Aim2EditFieldLabel.Position = [14 710 43 22];
            app.Aim2EditFieldLabel.Text = 'Ai(m^2)';

            % Create Aim2EditField
            app.Aim2EditField = uieditfield(app.RightPanel, 'numeric');
            app.Aim2EditField.Position = [103 710 85 22];

            % Create dcmEditFieldLabel
            app.dcmEditFieldLabel = uilabel(app.RightPanel);
            app.dcmEditFieldLabel.HorizontalAlignment = 'right';
            app.dcmEditFieldLabel.Position = [17 797 31 22];
            app.dcmEditFieldLabel.Text = 'd(cm)';

            % Create dcmEditField
            app.dcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.dcmEditField.Position = [103 797 84 22];

            % Create acmEditFieldLabel
            app.acmEditFieldLabel = uilabel(app.RightPanel);
            app.acmEditFieldLabel.HorizontalAlignment = 'right';
            app.acmEditFieldLabel.Position = [17 738 34 22];
            app.acmEditFieldLabel.Text = 'a(cm)';

            % Create acmEditField
            app.acmEditField = uieditfield(app.RightPanel, 'numeric');
            app.acmEditField.Position = [103 738 84 22];

            % Create hyokecmEditFieldLabel
            app.hyokecmEditFieldLabel = uilabel(app.RightPanel);
            app.hyokecmEditFieldLabel.HorizontalAlignment = 'right';
            app.hyokecmEditFieldLabel.Position = [410 770 63 22];
            app.hyokecmEditFieldLabel.Text = 'h yoke(cm)';

            % Create hyokecmEditField
            app.hyokecmEditField = uieditfield(app.RightPanel, 'numeric');
            app.hyokecmEditField.Position = [496 765 84 22];

            % Create WwcmEditFieldLabel
            app.WwcmEditFieldLabel = uilabel(app.RightPanel);
            app.WwcmEditFieldLabel.HorizontalAlignment = 'right';
            app.WwcmEditFieldLabel.Position = [216 766 46 22];
            app.WwcmEditFieldLabel.Text = 'Ww(cm)';

            % Create WwcmEditField
            app.WwcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.WwcmEditField.Position = [302 765 84 22];

            % Create HwcmEditFieldLabel
            app.HwcmEditFieldLabel = uilabel(app.RightPanel);
            app.HwcmEditFieldLabel.HorizontalAlignment = 'right';
            app.HwcmEditFieldLabel.Position = [212 741 47 22];
            app.HwcmEditFieldLabel.Text = 'Hw(cm)';

            % Create HwcmEditField
            app.HwcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.HwcmEditField.Position = [302 739 84 22];

            % Create DcmEditFieldLabel
            app.DcmEditFieldLabel = uilabel(app.RightPanel);
            app.DcmEditFieldLabel.HorizontalAlignment = 'right';
            app.DcmEditFieldLabel.Position = [14 765 35 22];
            app.DcmEditFieldLabel.Text = 'D(cm)';

            % Create DcmEditField
            app.DcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.DcmEditField.Position = [103 767 84 22];

            % Create CoreheighcmtEditFieldLabel
            app.CoreheighcmtEditFieldLabel = uilabel(app.RightPanel);
            app.CoreheighcmtEditFieldLabel.HorizontalAlignment = 'right';
            app.CoreheighcmtEditFieldLabel.Position = [389 738 91 22];
            app.CoreheighcmtEditFieldLabel.Text = 'Core heigh(cm)t';

            % Create CoreheighcmtEditField
            app.CoreheighcmtEditField = uieditfield(app.RightPanel, 'numeric');
            app.CoreheighcmtEditField.Position = [494 738 85 22];

            % Create CoreWidthcmEditFieldLabel
            app.CoreWidthcmEditFieldLabel = uilabel(app.RightPanel);
            app.CoreWidthcmEditFieldLabel.HorizontalAlignment = 'right';
            app.CoreWidthcmEditFieldLabel.Position = [395 711 86 22];
            app.CoreWidthcmEditFieldLabel.Text = 'Core Width(cm)';

            % Create CoreWidthcmEditField
            app.CoreWidthcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.CoreWidthcmEditField.Position = [495 710 85 22];

            % Create Awm2EditFieldLabel
            app.Awm2EditFieldLabel = uilabel(app.RightPanel);
            app.Awm2EditFieldLabel.HorizontalAlignment = 'right';
            app.Awm2EditFieldLabel.Position = [216 797 49 22];
            app.Awm2EditFieldLabel.Text = 'Aw(m^2)';

            % Create Awm2EditField
            app.Awm2EditField = uieditfield(app.RightPanel, 'numeric');
            app.Awm2EditField.Position = [302 797 84 22];

            % Create CoolingMethodEditFieldLabel
            app.CoolingMethodEditFieldLabel = uilabel(app.RightPanel);
            app.CoolingMethodEditFieldLabel.HorizontalAlignment = 'right';
            app.CoolingMethodEditFieldLabel.Position = [200 678 88 22];
            app.CoolingMethodEditFieldLabel.Text = 'Cooling Method';

            % Create CoolingMethodEditField
            app.CoolingMethodEditField = uieditfield(app.RightPanel, 'text');
            app.CoolingMethodEditField.Position = [301 677 85 22];

            % Create HVturnsEditFieldLabel
            app.HVturnsEditFieldLabel = uilabel(app.RightPanel);
            app.HVturnsEditFieldLabel.HorizontalAlignment = 'right';
            app.HVturnsEditFieldLabel.Position = [17 578 55 22];
            app.HVturnsEditFieldLabel.Text = 'H-V turns';

            % Create HVturnsEditField
            app.HVturnsEditField = uieditfield(app.RightPanel, 'numeric');
            app.HVturnsEditField.Position = [149 578 115 22];

            % Create LVturnsEditFieldLabel
            app.LVturnsEditFieldLabel = uilabel(app.RightPanel);
            app.LVturnsEditFieldLabel.HorizontalAlignment = 'right';
            app.LVturnsEditFieldLabel.Position = [17 459 49 22];
            app.LVturnsEditFieldLabel.Text = 'L-V turns';

            % Create LVturnsEditField
            app.LVturnsEditField = uieditfield(app.RightPanel, 'numeric');
            app.LVturnsEditField.Position = [149 459 114 22];

            % Create HVcurrentAEditFieldLabel
            app.HVcurrentAEditFieldLabel = uilabel(app.RightPanel);
            app.HVcurrentAEditFieldLabel.HorizontalAlignment = 'right';
            app.HVcurrentAEditFieldLabel.Position = [17 546 80 22];
            app.HVcurrentAEditFieldLabel.Text = 'H-V current(A)';

            % Create HVcurrentAEditField
            app.HVcurrentAEditField = uieditfield(app.RightPanel, 'numeric');
            app.HVcurrentAEditField.Position = [149 546 115 22];

            % Create LVcurrentAEditFieldLabel
            app.LVcurrentAEditFieldLabel = uilabel(app.RightPanel);
            app.LVcurrentAEditFieldLabel.HorizontalAlignment = 'right';
            app.LVcurrentAEditFieldLabel.Position = [17 433 76 22];
            app.LVcurrentAEditFieldLabel.Text = 'L-V current(A)';

            % Create LVcurrentAEditField
            app.LVcurrentAEditField = uieditfield(app.RightPanel, 'numeric');
            app.LVcurrentAEditField.Position = [149 433 114 22];

            % Create HVconductorcsamm2EditFieldLabel
            app.HVconductorcsamm2EditFieldLabel = uilabel(app.RightPanel);
            app.HVconductorcsamm2EditFieldLabel.HorizontalAlignment = 'right';
            app.HVconductorcsamm2EditFieldLabel.Position = [283 513 139 22];
            app.HVconductorcsamm2EditFieldLabel.Text = 'H-V conductor csa(mm2)';

            % Create HVconductorcsamm2EditField
            app.HVconductorcsamm2EditField = uieditfield(app.RightPanel, 'numeric');
            app.HVconductorcsamm2EditField.Position = [466 513 116 22];

            % Create LVconductorcsamm2EditFieldLabel
            app.LVconductorcsamm2EditFieldLabel = uilabel(app.RightPanel);
            app.LVconductorcsamm2EditFieldLabel.HorizontalAlignment = 'right';
            app.LVconductorcsamm2EditFieldLabel.Position = [283 404 137 22];
            app.LVconductorcsamm2EditFieldLabel.Text = 'L-V conductor csa(mm2)';

            % Create LVconductorcsamm2EditField
            app.LVconductorcsamm2EditField = uieditfield(app.RightPanel, 'numeric');
            app.LVconductorcsamm2EditField.Position = [466 407 116 22];

            % Create LVinnerdiametercmEditFieldLabel
            app.LVinnerdiametercmEditFieldLabel = uilabel(app.RightPanel);
            app.LVinnerdiametercmEditFieldLabel.HorizontalAlignment = 'right';
            app.LVinnerdiametercmEditFieldLabel.Position = [287 459 124 22];
            app.LVinnerdiametercmEditFieldLabel.Text = 'L-V inner diameter(cm)';

            % Create LVinnerdiametercmEditField
            app.LVinnerdiametercmEditField = uieditfield(app.RightPanel, 'numeric');
            app.LVinnerdiametercmEditField.Position = [466 459 116 22];

            % Create XeqEditFieldLabel
            app.XeqEditFieldLabel = uilabel(app.RightPanel);
            app.XeqEditFieldLabel.HorizontalAlignment = 'right';
            app.XeqEditFieldLabel.Position = [14 259 24 22];
            app.XeqEditFieldLabel.Text = 'Xeq';

            % Create XeqEditField
            app.XeqEditField = uieditfield(app.RightPanel, 'numeric');
            app.XeqEditField.Position = [103 259 85 22];

            % Create LVHeightcmEditFieldLabel
            app.LVHeightcmEditFieldLabel = uilabel(app.RightPanel);
            app.LVHeightcmEditFieldLabel.HorizontalAlignment = 'right';
            app.LVHeightcmEditFieldLabel.Position = [16 406 82 22];
            app.LVHeightcmEditFieldLabel.Text = 'L-V Height(cm)';

            % Create LVHeightcmEditField
            app.LVHeightcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.LVHeightcmEditField.Position = [149 406 114 22];

            % Create FormarHeightcmEditFieldLabel
            app.FormarHeightcmEditFieldLabel = uilabel(app.RightPanel);
            app.FormarHeightcmEditFieldLabel.HorizontalAlignment = 'right';
            app.FormarHeightcmEditFieldLabel.Position = [14 364 105 22];
            app.FormarHeightcmEditFieldLabel.Text = 'Formar Height (cm)';

            % Create FormarHeightcmEditField
            app.FormarHeightcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.FormarHeightcmEditField.Position = [149 364 115 22];

            % Create FormarwidthcmEditFieldLabel
            app.FormarwidthcmEditFieldLabel = uilabel(app.RightPanel);
            app.FormarwidthcmEditFieldLabel.HorizontalAlignment = 'right';
            app.FormarwidthcmEditFieldLabel.Position = [286 364 96 22];
            app.FormarwidthcmEditFieldLabel.Text = 'Formar width(cm)';

            % Create FormarwidthcmEditField
            app.FormarwidthcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.FormarwidthcmEditField.Position = [466 364 116 22];

            % Create HVHeightcmEditFieldLabel
            app.HVHeightcmEditFieldLabel = uilabel(app.RightPanel);
            app.HVHeightcmEditFieldLabel.HorizontalAlignment = 'right';
            app.HVHeightcmEditFieldLabel.Position = [287 578 83 22];
            app.HVHeightcmEditFieldLabel.Text = 'H-V Height(cm)';

            % Create HVHeightcmEditField
            app.HVHeightcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.HVHeightcmEditField.Position = [466 578 114 22];

            % Create HVWidthcmEditFieldLabel
            app.HVWidthcmEditFieldLabel = uilabel(app.RightPanel);
            app.HVWidthcmEditFieldLabel.HorizontalAlignment = 'right';
            app.HVWidthcmEditFieldLabel.Position = [287 543 80 22];
            app.HVWidthcmEditFieldLabel.Text = 'H-V Width(cm)';

            % Create HVWidthcmEditField
            app.HVWidthcmEditField = uieditfield(app.RightPanel, 'numeric');
            app.HVWidthcmEditField.Position = [466 543 114 22];

            % Create ZbaseEditFieldLabel
            app.ZbaseEditFieldLabel = uilabel(app.RightPanel);
            app.ZbaseEditFieldLabel.HorizontalAlignment = 'right';
            app.ZbaseEditFieldLabel.Position = [448 290 47 22];
            app.ZbaseEditFieldLabel.Text = 'Zbase';

            % Create ZbaseEditField
            app.ZbaseEditField = uieditfield(app.RightPanel, 'numeric');
            app.ZbaseEditField.Position = [505 290 71 22];

            % Create HVconductortypeEditField_2Label
            app.HVconductortypeEditField_2Label = uilabel(app.RightPanel);
            app.HVconductortypeEditField_2Label.HorizontalAlignment = 'right';
            app.HVconductortypeEditField_2Label.Position = [17 513 107 22];
            app.HVconductortypeEditField_2Label.Text = 'H-V conductor type';

            % Create HVconductortypeEditField_2
            app.HVconductortypeEditField_2 = uieditfield(app.RightPanel, 'text');
            app.HVconductortypeEditField_2.Position = [149 513 115 22];

            % Create LVconductortypeEditFieldLabel
            app.LVconductortypeEditFieldLabel = uilabel(app.RightPanel);
            app.LVconductortypeEditFieldLabel.HorizontalAlignment = 'right';
            app.LVconductortypeEditFieldLabel.Position = [287 433 103 22];
            app.LVconductortypeEditFieldLabel.Text = 'L-V conductor type';

            % Create LVconductortypeEditField
            app.LVconductortypeEditField = uieditfield(app.RightPanel, 'text');
            app.LVconductortypeEditField.Position = [466 433 116 22];

            % Create IcoreAEditFieldLabel
            app.IcoreAEditFieldLabel = uilabel(app.RightPanel);
            app.IcoreAEditFieldLabel.HorizontalAlignment = 'right';
            app.IcoreAEditFieldLabel.Position = [17 229 47 22];
            app.IcoreAEditFieldLabel.Text = 'I core(A)';

            % Create IcoreAEditField
            app.IcoreAEditField = uieditfield(app.RightPanel, 'numeric');
            app.IcoreAEditField.Position = [149 229 115 22];

            % Create InoloadAEditFieldLabel
            app.InoloadAEditFieldLabel = uilabel(app.RightPanel);
            app.InoloadAEditFieldLabel.HorizontalAlignment = 'right';
            app.InoloadAEditFieldLabel.Position = [394 259 67 22];
            app.InoloadAEditFieldLabel.Text = 'I no load(A)';

            % Create InoloadAEditField
            app.InoloadAEditField = uieditfield(app.RightPanel, 'numeric');
            app.InoloadAEditField.Position = [505 259 74 22];

            % Create RcoreohmEditFieldLabel
            app.RcoreohmEditFieldLabel = uilabel(app.RightPanel);
            app.RcoreohmEditFieldLabel.HorizontalAlignment = 'right';
            app.RcoreohmEditFieldLabel.Position = [287 229 68 22];
            app.RcoreohmEditFieldLabel.Text = 'R core(ohm)';

            % Create RcoreohmEditField
            app.RcoreohmEditField = uieditfield(app.RightPanel, 'numeric');
            app.RcoreohmEditField.Position = [466 229 113 22];

            % Create ImagtAEditFieldLabel
            app.ImagtAEditFieldLabel = uilabel(app.RightPanel);
            app.ImagtAEditFieldLabel.HorizontalAlignment = 'right';
            app.ImagtAEditFieldLabel.Position = [212 259 71 22];
            app.ImagtAEditFieldLabel.Text = 'I magt(A)';

            % Create ImagtAEditField
            app.ImagtAEditField = uieditfield(app.RightPanel, 'numeric');
            app.ImagtAEditField.Position = [302 259 84 22];

            % Create CopperlosswattEditFieldLabel
            app.CopperlosswattEditFieldLabel = uilabel(app.RightPanel);
            app.CopperlosswattEditFieldLabel.HorizontalAlignment = 'right';
            app.CopperlosswattEditFieldLabel.Position = [287 194 95 22];
            app.CopperlosswattEditFieldLabel.Text = 'Copper loss(watt)';

            % Create CopperlosswattEditField
            app.CopperlosswattEditField = uieditfield(app.RightPanel, 'numeric');
            app.CopperlosswattEditField.Position = [466 194 113 22];

            % Create ExpectedefficiencyEditFieldLabel
            app.ExpectedefficiencyEditFieldLabel = uilabel(app.RightPanel);
            app.ExpectedefficiencyEditFieldLabel.HorizontalAlignment = 'right';
            app.ExpectedefficiencyEditFieldLabel.Position = [14 194 120 22];
            app.ExpectedefficiencyEditFieldLabel.Text = 'Expected efficiency %';

            % Create ExpectedefficiencyEditField
            app.ExpectedefficiencyEditField = uieditfield(app.RightPanel, 'numeric');
            app.ExpectedefficiencyEditField.Position = [149 194 115 22];

            % Create TankMeasurementsaDimensionsLabel
            app.TankMeasurementsaDimensionsLabel = uilabel(app.RightPanel);
            app.TankMeasurementsaDimensionsLabel.FontWeight = 'bold';
            app.TankMeasurementsaDimensionsLabel.FontColor = [0.4941 0.1843 0.5569];
            app.TankMeasurementsaDimensionsLabel.Position = [32 144 534 22];
            app.TankMeasurementsaDimensionsLabel.Text = '.............................................. Tank Measurementsa & Dimensions ..............................................';

            % Create TanklengthXmEditFieldLabel
            app.TanklengthXmEditFieldLabel = uilabel(app.RightPanel);
            app.TanklengthXmEditFieldLabel.HorizontalAlignment = 'right';
            app.TanklengthXmEditFieldLabel.Position = [-2 114 97 22];
            app.TanklengthXmEditFieldLabel.Text = 'Tank length''X''(m)';

            % Create TanklengthXmEditField
            app.TanklengthXmEditField = uieditfield(app.RightPanel, 'numeric');
            app.TanklengthXmEditField.Position = [103 114 85 22];

            % Create TankwidthYmEditFieldLabel
            app.TankwidthYmEditFieldLabel = uilabel(app.RightPanel);
            app.TankwidthYmEditFieldLabel.HorizontalAlignment = 'right';
            app.TankwidthYmEditFieldLabel.Position = [193 114 93 22];
            app.TankwidthYmEditFieldLabel.Text = 'Tank width''Y''(m)';

            % Create TankwidthYmEditField
            app.TankwidthYmEditField = uieditfield(app.RightPanel, 'numeric');
            app.TankwidthYmEditField.Position = [302 114 84 22];

            % Create TankheightZmEditFieldLabel
            app.TankheightZmEditFieldLabel = uilabel(app.RightPanel);
            app.TankheightZmEditFieldLabel.HorizontalAlignment = 'right';
            app.TankheightZmEditFieldLabel.Position = [384 114 97 22];
            app.TankheightZmEditFieldLabel.Text = 'Tank height''Z''(m)';

            % Create TankheightZmEditField
            app.TankheightZmEditField = uieditfield(app.RightPanel, 'numeric');
            app.TankheightZmEditField.Position = [505 114 76 22];

            % Create Tankdissipationareainm2EditFieldLabel
            app.Tankdissipationareainm2EditFieldLabel = uilabel(app.RightPanel);
            app.Tankdissipationareainm2EditFieldLabel.HorizontalAlignment = 'right';
            app.Tankdissipationareainm2EditFieldLabel.Position = [287 82 162 22];
            app.Tankdissipationareainm2EditFieldLabel.Text = 'Tank dissipation area in m^2%';

            % Create Tankdissipationareainm2EditField
            app.Tankdissipationareainm2EditField = uieditfield(app.RightPanel, 'numeric');
            app.Tankdissipationareainm2EditField.Position = [466 82 116 22];

            % Create MaximumtemperatureriseCLabel
            app.MaximumtemperatureriseCLabel = uilabel(app.RightPanel);
            app.MaximumtemperatureriseCLabel.HorizontalAlignment = 'right';
            app.MaximumtemperatureriseCLabel.Position = [14 82 169 22];
            app.MaximumtemperatureriseCLabel.Text = 'Maximum temperature rise (°C)';

            % Create MaximumtemperatureriseCEditField
            app.MaximumtemperatureriseCEditField = uieditfield(app.RightPanel, 'numeric');
            app.MaximumtemperatureriseCEditField.Position = [194 82 69 22];

            % Create CoolingtubenecessityEditFieldLabel
            app.CoolingtubenecessityEditFieldLabel = uilabel(app.RightPanel);
            app.CoolingtubenecessityEditFieldLabel.HorizontalAlignment = 'right';
            app.CoolingtubenecessityEditFieldLabel.Position = [1 47 126 22];
            app.CoolingtubenecessityEditFieldLabel.Text = 'Cooling tube necessity';

            % Create CoolingtubenecessityEditField
            app.CoolingtubenecessityEditField = uieditfield(app.RightPanel, 'text');
            app.CoolingtubenecessityEditField.Position = [133 47 195 22];

            % Create DesignFeasibilityEditFieldLabel
            app.DesignFeasibilityEditFieldLabel = uilabel(app.RightPanel);
            app.DesignFeasibilityEditFieldLabel.HorizontalAlignment = 'right';
            app.DesignFeasibilityEditFieldLabel.Position = [331 47 100 22];
            app.DesignFeasibilityEditFieldLabel.Text = 'Design Feasibility';

            % Create DesignFeasibilityEditField
            app.DesignFeasibilityEditField = uieditfield(app.RightPanel, 'text');
            app.DesignFeasibilityEditField.Position = [437 47 145 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Transformer_Simulator_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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