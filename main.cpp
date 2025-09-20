/*
 * This file is part of the ProVANT simulator project.
 * Licensed under the terms of the MIT open source license. More details at
 * https://github.com/Guiraffo/ProVANT-Simulator/blob/master/LICENSE.md
 */
/**
 * @file This file contains the implementation of the MpcQuadtilt class.
 *
 * @author Insert Your Name Here
 */

#include <control_strategies_base/icontroller.hpp>
#include <iostream>
#include <Eigen/Eigen>
#include "simulator_msgs/Sensor.h"
#include <ros/ros.h>
#include <vector>
#include "modelo.hpp"
#include <std_srvs/Empty.h>
#include <cmath>
using namespace Eigen;

static const double amost = 0.01; // amostragem
static const int N = 60; // Horizonte de Predição

// Parâmetros simbólicos
static const int estados = 16;
static const int entradas = 4;

// Restrições e limites
std::vector<double> lbg(estados*(N+1), 0);
std::vector<double> ubg(estados*(N+1), 0);

std::vector<double> lbx(estados*(N+1) + entradas*N);
std::vector<double> ubx(estados*(N+1) + entradas*N);

//DEFINIR ESTADOS x=[e, v, q, ohm]' \in R^13
SX pose = SX::sym("pose",8,1);
SX heligiro = SX::sym("heligiro",8,1);
std::vector<SX> states = {pose, heligiro};


//DIFINIR ENTRADAS u=[f1, f2, f3, f4] \in R^4
SX f1 = SX::sym("f1");
SX f2 = SX::sym("f2");
SX f3 = SX::sym("f3");
SX f4 = SX::sym("f4");
std::vector<SX> inputs_func = {f1,f2,f3,f4};


//Vectores dos estados a ser calculados
SX U = SX::sym("U", entradas, N);
SX X = SX::sym("X", estados, N+1); //toods os estados ao longo do N
SX X0 = SX::sym("X0", estados); //estado inicial
SX X_r = SX::sym("X_r", estados); //estado de referencia
//SX P = SX::sym("P", estados+estados);
int num_err = 0;

// Funcional de custo
SX J = 0;

// Criação da matriz Q e Qh usando CasADi
SX Q = SX::zeros(8,8);
SX Qh = SX::zeros(8,8);

// Cria uma matriz 4x4 com zeros
SX R = SX::zeros(4,4);

static const int K_psi = 220;
static const int zeta = 1.045;

//Inicializa restrições
SX g;

//cria solver
Function solver;


// Parâmetros iniciais
DM x0;
DM x_0;


//Referencia
DM x_r;


//Define entrada de controle inicial
DM u_0;
static const int NINPUTS = 4;
static const int NSTATES = 16;

//Cria variavel usada para guardar as infromações do solver
std::map<std::string, DM> args;


ros::NodeHandle nh;

    // Clientes de serviços para pausar e despausar a física do Gazebo
ros::ServiceClient pause_physics_client = nh.serviceClient<std_srvs::Empty>("/gazebo/pause_physics");
ros::ServiceClient unpause_physics_client = nh.serviceClient<std_srvs::Empty>("/gazebo/unpause_physics");
ros::ServiceClient step_physics_client = nh.serviceClient<std_srvs::Empty>("/gazebo/step");

    // Declarar o serviço uma única vez fora do loop
 std_srvs::Empty srv;
//entrada de controle



class MpcQuadtilt : public Icontroller
{
private:
Eigen::VectorXd Input;
Eigen::VectorXd X_STATES;
Eigen::VectorXd Xref;
Eigen::VectorXd Erro;
Eigen::VectorXd Erro_ref;
Eigen::VectorXd Erro_state;
std::vector<double> alvo;           // Alvo como membro da classe
std::vector<double> xr_init;        // Vetor de referência
Eigen::Quaterniond q_d;
static bool has_landed;             // Variável estática
static ros::Time landed_time; // Armazena o tempo em que o drone pousou

public: MpcQuadtilt() :  Xref(NSTATES),X_STATES(NSTATES), Erro(NSTATES), Input(NINPUTS)
        {
  config();
        }



    virtual ~MpcQuadtilt() {}




void pauseGazebo(ros::ServiceClient& pause_physics_client, std_srvs::Empty& srv) {

   pause_physics_client.call(srv);

}

void unpauseGazebo(ros::ServiceClient& unpause_physics_client, std_srvs::Empty& srv) {
   unpause_physics_client.call(srv);
}

void stepGazebo(ros::ServiceClient& step_physics_client, std_srvs::Empty& srv) {
    step_physics_client.call(srv);
}
    void conversao()
    {

        //X = [pose, heligiro]
        // Limites para a pose primária
        x0(0) = X_STATES(0);
        x0(1) = X_STATES(1);
        x0(2) = X_STATES(2);
        x0(3) = X_STATES(3);

        // Limites para a pose dual
        x0(4) = X_STATES(4);
        x0(5) = X_STATES(5);
        x0(6) = X_STATES(6);
        x0(7) = X_STATES(7);

        // Limites para o heligiro primário
        x0(8) = X_STATES(8);
        x0(9) = X_STATES(9);
        x0(10) = X_STATES(10);
        x0(11) = X_STATES(11);

        // Limites para o heligiro dual
        x0(12) = X_STATES(12);
        x0(13) = X_STATES(13);
        x0(14) = X_STATES(14);
        x0(15) = X_STATES(15);

        x0 =  x0;
        std::cout << "x0" << x0 << std::endl;

    }
    void bounds()
    {

        for (int i = 0; i < N+1; i++)
        {
            // Limites para a pose primária
                   lbx[i * estados + 0] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 0] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 1] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 1] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 2] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 2] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 3] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 3] = std::numeric_limits<double>::infinity();
                   
                   // Limites para a pose dual
                   lbx[i * estados + 4] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 4] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 5] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 5] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 6] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 6] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 7] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 7] = std::numeric_limits<double>::infinity();
                   
                   // Limites para o heligiro primário
                   lbx[i * estados + 8] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 8] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 9] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 9] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 10] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 10] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 11] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 11] = std::numeric_limits<double>::infinity();

                   // Limites para o heligiro dual
                   lbx[i * estados + 12] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 12] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 13] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 13] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 14] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 14] = std::numeric_limits<double>::infinity();

                   lbx[i * estados + 15] = -std::numeric_limits<double>::infinity();
                   ubx[i * estados + 15] = std::numeric_limits<double>::infinity();



        }
        for (int i = 0; i < N; i++)
        {
            lbx[estados * (N + 1) + i * entradas + 0] = 0;       // f1 lower bound
                    ubx[estados * (N + 1) + i * entradas + 0] = 12.3;    // f1 upper bound

                    lbx[estados * (N + 1) + i * entradas + 1] = 0;       // f2 lower bound
                    ubx[estados * (N + 1) + i * entradas + 1] = 12.3;    // f2 upper bound

                    lbx[estados * (N + 1) + i * entradas + 2] = 0;       // f3 lower bound
                    ubx[estados * (N + 1) + i * entradas + 2] = 12.3;    // f3 upper bound

                    lbx[estados * (N + 1) + i * entradas + 3] = 0;       // f4 lower bound
                    ubx[estados * (N + 1) + i * entradas + 3] = 12.3;    // f4 upper bound
        }

    }

    void tunner()
    {


  // Modifica os elementos da diagonal da matriz Q com os valores fornecidos
    Q(0, 0) = 0;
    Q(1, 1) = 50;
    Q(2, 2) = 50;
    Q(3, 3) = 360;
    Q(4, 4) = 0;
    Q(5, 5) = 680;
    Q(6, 6) = 680;
    Q(7, 7) = 490;

  // Modifica os elementos da diagonal da matriz Qh com os valores fornecidos
    Qh(0, 0) = 0;
    Qh(1, 1) = 10;
    Qh(2, 2) = 10;
    Qh(3, 3) = 10;
    Qh(4, 4) = 0;
    Qh(5, 5) = 65;
    Qh(6, 6) = 65;
    Qh(7, 7) = 45;

    // Modifica os elementos da diagonal da matriz R com os valores fornecidos
    R(0, 0) = 0.05; // f1
    R(1, 1) = 0.05; // f2
    R(2, 2) = 0.05; // f3
    R(3, 3) = 0.05; // f4

    }



    void config()
    {

        tunner();
        g = X(Slice(), 0) - X0;
        SX p_dot = SX::zeros(8, 1);
        SX aak = SX::zeros(8, 1);
        SX wb = SX::zeros(8, 1);
        SX x_outro = SX::zeros(8, 1);

        // SX q1= SX::vertcat({1, 2, 3, 4, 5, 6, 7, 8});
        // SX q2 = SX::vertcat({8, 7, 6, 5, 4, 3, 2, 1});
        // SX unit = SX::vertcat({0.999494449756450, 0.0307811844561152, -0.00532928744955115, 0.00591289140851692, -0.0506967755136625, 1.55958164242684, 0.552845955329595, 0.949053544215771});

        // SX test = normalize_qd_casadi(unit);
        // std::cout << "normalize: " << test << std::endl; //pode excluir


    for (int k = 0; k < N; k++)
    {

        SX X_k = X(Slice(), k);
        SX U_k = U(Slice(), k);
        std::cout << k << std::endl;
        


        SX P_1=vertcat_with_check({X0, X_r}, {"X0", "X_r"});

        J += compute_J_qd(Q, Qh, R, X0, X_r, X_k, U_k);
        // std::cout << "J: " << J << std::endl; //pode excluir

        SX st_next = X(Slice(), k+1);  // Next state
        SX st_next_qd;
        std::tie(st_next_qd, p_dot, aak, wb, x_outro) = compute_next_qd(X_k, U_k, p_dot, aak, wb, x_outro, amost);
        // std::cout << "st_next_qd: " << st_next_qd << std::endl; //pode excluir

        g = vertcat_with_check({g, st_next - st_next_qd}, {"g", "st_next - st_next_qd"});
        // std::cout << "g: " << g << std::endl; //pode excluir

    }
    g = simplify(g);


        // Variáveis de otimização
        SX opt_variables = vertcat(reshape(X, estados*(N+1), 1), reshape(U, entradas*N, 1));
        SX P = vertcat(X0, X_r);


        //std::cout << "Dimensão de X: " << X.size1() << " x " << X.size2() << std::endl; // Deve ser 13 x 51
        //std::cout << "Dimensão de U: " << U.size1() << " x " << U.size2() << std::endl; // Deve ser 4 x 50

        // Criação do problema NLP
        std::map<std::string, SX> nlp_prob =
        {
            {"f", J},
            {"x", opt_variables},
            {"g", g},
            {"p", P}
        };


        std::string solver_name = "ipopt";
        Dict nlp_opts;
        nlp_opts["expand"] = true;
        nlp_opts["ipopt.max_iter"] = 2000;
        nlp_opts["ipopt.print_level"] = 0;
        nlp_opts["print_time"] = 0;
        nlp_opts["ipopt.acceptable_tol"] =  1e-8;
        nlp_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;


        // Criação do solver
        solver = nlpsol("solver", "ipopt", nlp_prob, nlp_opts);
        bounds();



        // Estado inicial
        SX pos_init = SX::vertcat({0, 0, 0, 0});
        SX ori_init = SX::vertcat({1, 0, 0, 0});
        DM vel_init = DM::vertcat({0, 0, 0, 0});
        DM ome_init = DM::vertcat({0, 0.1, 0.1, 0.1});

        DM dual_pose_init = 0.5 * DM(produto_q(pos_init, ori_init));
        DM pos_completa_init = DM::vertcat({DM(ori_init), dual_pose_init});
        DM heligiro_init = DM::vertcat({vel_init, ome_init});

        DM x0_init = DM::vertcat({pos_completa_init, heligiro_init});
        // Chamar a função Init
        x0 = x0_init;
        //x0 = Init(x0_init);
        x_0 = repmat(x0, 1, N+1);
        DM Ueq = DM::vertcat({0.0, 0.0, 0.0,0.0});
        u_0 = repmat(Ueq, 1, N);



        // Modificar valores específicos conforme necessário
        //alvo = {3,3.8,1.2};
        alvo = {1,1,1};

        double delta_x = alvo[0] - static_cast<double>(x0(0).scalar());
        double delta_y = alvo[1] - static_cast<double>(x0(1).scalar());
        double yaw_mexe = std::atan2(delta_y, delta_x);

        std::cout << "yaw_mexe: " << yaw_mexe << std::endl;
        std::cout << "yaw_mexe (em graus): " << yaw_mexe * 180 / M_PI << std::endl;


        q_d = eul2quat(0, 0, yaw_mexe);

        std::cout << q_d.w() << " " << q_d.x() << " " << q_d.y() << " " << q_d.z() << std::endl;
        SX p_des = SX::vertcat({SX(0), SX(alvo[0]), SX(alvo[1]), SX(alvo[2])});
        SX q_des = SX::vertcat({SX(q_d.w()), SX(q_d.x()), SX(q_d.y()), SX(q_d.z())});
        DM dual_pose_des = 0.5 * DM(produto_q(p_des, q_des));

    xr_init = {
    q_d.w(), q_d.x(), q_d.y(), q_d.z(),     // Orientação desejada
    static_cast<double>(dual_pose_des(0)),
    static_cast<double>(dual_pose_des(1)),
    static_cast<double>(dual_pose_des(2)),
    static_cast<double>(dual_pose_des(3)),
    0.0, 0.0, 0.0, 0.0,                     // Velocidade angular desejada
    0.0, 0.0, 0.0, 0.0                      // Velocidade linear desejada
};

        //DM xr_init = DM{4.4, 5.4, 2.1, 0, 0, 0, 0.9421, 0, 0, 0.3352, 0, 0, 0};

        // Chamar a função Init
        x_r = xr_init;


        args["p"] = vertcat(x0, x_r);
        args["x0"] = vertcat(reshape(x_0, estados*(N+1), 1), reshape(u_0, entradas*N, 1));
        args["lbx"] = lbx;
        args["ubx"] = ubx;
        args["lbg"] = lbg;
        args["ubg"] = ubg;
        std::cout << "Espere2\n"<< std::endl;

    }

    SX prox_pos = SX::zeros(4, 1);
    SX prox_ori = SX::zeros(4, 1);
    SX prox_vel = SX::zeros(4, 1);
    SX prox_ome = SX::zeros(4, 1);
    SX prox_dual_pose = SX::zeros(4, 1);

    std::vector<double> execute(simulator_msgs::SensorArray arraymsg)
    {



        pauseGazebo(pause_physics_client, srv);
        simulator_msgs::Sensor msg;
        bool found = false;
        for (int i = 0; i < arraymsg.values.size(); i++)
        {
            if (arraymsg.values.at(i).name == "Estados")
            {
                msg = arraymsg.values.at(i);
                found = true;
                std::cout << "Estados obtidos:\n" << std::endl;
                std::cout << arraymsg.values << std::endl;

                break;
            }
        }
        if (!found)
        {
            // In case of error, report the problem in a ROS_LOG, and returns an
            // empty array
            ROS_FATAL("[lqr_vant5] State vector not found.");
            std::vector<double> out(Input.data(), Input.data() + Input.size());
            std::cout << "erro ao receber msgs"<< std::endl;
            return out;
        }

 // Get the current simulation time and calculate the reference trajectory.
    ros::Time time = ros::Time::now();
    double tempo = time.toSec();

    // Atualizar o X_STATES e o Xref para atualizar os vetores para a proxima iteracaos
    Xref << xr_init[0], xr_init[1], xr_init[2], xr_init[3],
             xr_init[4], xr_init[5], xr_init[6], xr_init[7],
             xr_init[8], xr_init[9], xr_init[10], xr_init[11],
             xr_init[12], xr_init[13], xr_init[14], xr_init[15];
    
    // Read the values of the state vector from the message 
    
    prox_pos = SX::vertcat({0, msg.values.at(0), msg.values.at(1), msg.values.at(2)});
    prox_ori = SX::vertcat({msg.values.at(21), msg.values.at(18), msg.values.at(19), msg.values.at(20)});
    prox_vel = SX::vertcat({0, msg.values.at(6), msg.values.at(7), msg.values.at(8)});
    prox_ome = SX::vertcat({0, msg.values.at(15), msg.values.at(16), msg.values.at(17)});

    prox_dual_pose = 0.5 * produto_q(prox_pos, prox_ori);
    // prox_pose = SX::vertcat({prox_ori, prox_dual_pose}); //pode excluir
    prox_ome = ad_qd_casadi(conj_qd_casadi(SX::vertcat({prox_ori, SX::zeros(4,1)})), SX::vertcat({prox_ome, SX::zeros(4,1)}));
    // prox_heligiro = SX::vertcat({prox_ome, prox_vel}); //pode excluir

    // std::cout << "prox_pos: " << prox_pos << std::endl; //pode excluir
    // std::cout << "prox_ori: " << prox_ori << std::endl; //pode excluir
    // std::cout << "prox_dual_pose: " << prox_dual_pose << std::endl; //pode excluir
    // std::cout << "prox_ome: " << prox_ome << std::endl;//pode excluir
    // std::cout << "prox_vel: " << prox_vel << std::endl;//pode excluir

    // std::cout << "X_STATES antes" << X_STATES << std::endl;//pode excluir

    X_STATES << static_cast<double>(prox_ori(0).scalar()), static_cast<double>(prox_ori(1).scalar()), static_cast<double>(prox_ori(2).scalar()), static_cast<double>(prox_ori(3).scalar()),
                static_cast<double>(prox_dual_pose(0).scalar()), static_cast<double>(prox_dual_pose(1).scalar()), static_cast<double>(prox_dual_pose(2).scalar()), static_cast<double>(prox_dual_pose(3).scalar()),
                static_cast<double>(prox_ome(0).scalar()), static_cast<double>(prox_ome(1).scalar()), static_cast<double>(prox_ome(2).scalar()), static_cast<double>(prox_ome(3).scalar()),
                static_cast<double>(prox_vel(0).scalar()), static_cast<double>(prox_vel(1).scalar()), static_cast<double>(prox_vel(2).scalar()), static_cast<double>(prox_vel(3).scalar());
    
    std::cout << "X_STATES depois" << X_STATES << std::endl;

    // Calculate the error vector of the system
    // Erro = X_STATES - Xref;
    Erro_state << static_cast<double>(prox_ori(1).scalar()),
        static_cast<double>(prox_ori(2).scalar()),
        static_cast<double>(prox_ori(3).scalar());

    Erro_ref << alvo[0],
        alvo[1],
        alvo[2];

    Erro = Erro_state - Erro_ref;

    std::cout << "Erro: " << Erro << std::endl;


    // Calcula a norma apenas dos primeiros 3 elementos de Erro
    double erro_norma = Erro.head(10).norm();

    //double limite_erro = 0.07; // Defina o limite de erro para considerar que o alvo foi alcançado

    // Se o drone já pousou, mantém as entradas mínimas


        // Verifica se o drone está próximo do alvo para pousar
        double limite_erro = 0.07;

    std::cout << "Norma do erro: " << erro_norma << std::endl;


        conversao();
        // saida é x0=[x,y,z,vx,vy,vz,roll_dot,pitch_dot,yaw_dot,roll,pitch,yaw]
        x_0 = repmat(x0, 1, N+1);


        args["p"] = vertcat(x0, x_r);
        args["x0"] = vertcat(reshape(x_0, estados*(N+1), 1), reshape(u_0, entradas*N, 1));


        std::map<std::string, DM> sol = solver(args);

    std::map<std::string, casadi::GenericType> stats = solver.stats();
    bool success = static_cast<bool>(stats["success"]);
    if(success){
        //essa é a parte da simulação usando o while no matlab
        DM memoria = sol.at("x");
        u_0 = reshape(memoria(Slice(estados*(N+1), memoria.size1())), entradas, N);
        DM saida = u_0(Slice(),0);

        if (!has_landed)
        {
            Input(0) = double(saida(0));
            Input(1) = double(saida(1));
            Input(2) = double(saida(2));
            Input(3) = double(saida(3));

            double limite_erro = 0.07;
            if (erro_norma < limite_erro)
            {
                std::cout << "Drone chegou ao alvo. Simulando pouso." << std::endl;
                Input.setConstant(0.0); // Mantém entradas em 0 inicialmente
                has_landed = true;      // Marca como pousado
                landed_time = ros::Time::now(); // Registra o tempo do pouso
            }
        }
        else
        {
            // Calcula o tempo desde o pouso
            double elapsed_time = (ros::Time::now() - landed_time).toSec();

            if (elapsed_time < 0.5)
            {
                // Mantém as entradas em 0 por meio segundo
                Input.setConstant(0.0);
            }
            else
            {
                // Após meio segundo, pode definir outras entradas ou mantê-las constantes
                Input.setConstant(0.1); // Exemplo: manter entradas mínimas
            }
        }

        //std::cout << "Input analise" << Input << std::endl;
        num_err = 0;


     }else{
       num_err = num_err+1;
     std::cout << "Solver nao conseguiu resolver o problema de otimização:\n"<<num_err<<std::endl;
     }



    //u_0 = SX::vertcat({u_0(Slice(), Slice(1, N)), u_0(Slice(), Slice(N-1, N))}).T();
     u_0 = horzcat(u_0(Slice(), Slice(1, N)), u_0(Slice(), N-1));
    //u_0 = horzcat(u_0(Slice(), Slice(1, N)), u_0(Slice(), Slice(N-1, N)));

     std::cout << "u_0 valores" << u_0 << std::endl;

       unpauseGazebo(unpause_physics_client, srv);



        std::vector<double> out(Input.data(), Input.data() + Input.size());

        return out;
    }

    std::vector<double> Reference()
    {
        std::vector<double> out(Xref.data(), Xref.data() + Xref.rows() * Xref.cols());
        return out;
    }

    std::vector<double> Error()
    {
        std::vector<double> out(Erro.data(), Erro.data() + Erro.rows() * Erro.cols());
        return out;
    }

    std::vector<double> State()
    {
        std::vector<double> out(X_STATES.data(), X_STATES.data() + X_STATES.rows() * X_STATES.cols());
        return out;
    }
};

bool MpcQuadtilt::has_landed = false;
ros::Time MpcQuadtilt::landed_time = ros::Time(0);


extern "C"
{
    Icontroller *create(void)
    {
        return new MpcQuadtilt;
    }
    void destroy(Icontroller *p)
    {
        delete p;
    }
}
