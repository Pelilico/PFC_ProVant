#ifndef MODELO
#define MODELO
#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include <casadi/casadi.hpp>
#include <fstream>
#include <cstdlib>
#include <chrono>
using namespace std;
using namespace casadi;


// Parâmetros físicos globais
const double m = 2.24;       // Massa (kg)
const double d = 0.332;      // Comprimento do braço (m)
const double gr = 9.81;      // Aceleração devido à gravidade (m/s^2)
const double Ixx = 0.0384;   // Momento de inércia em torno do eixo x (kg·m^2)
const double Iyy = 0.0384;   // Momento de inércia em torno do eixo y (kg·m^2)
const double Izz = 0.0615;   // Momento de inércia em torno do eixo z (kg·m^2)
const double Ixxinv = 26.0417;
const double Iyyinv = 26.0417;
const double Izzinv = 16.2602;
const double ctau_f = 0.0179;  // Coeficiente de torque
const SX e_z = SX::vertcat({0, 0, 1}); // Vetor e_z
// Matriz de controle global
const SX mat_aux_con = []() {
    SX mat = SX::zeros(4, 4);
    mat(0, 0) = 1;   mat(0, 1) = 1;   mat(0, 2) = 1;   mat(0, 3) = 1;
    mat(1, 0) = 0;   mat(1, 1) = d;   mat(1, 2) = 0;   mat(1, 3) = -d;
    mat(2, 0) = -d;  mat(2, 1) = 0;   mat(2, 2) = d;   mat(2, 3) = 0;
    mat(3, 0) = ctau_f; mat(3, 1) = -ctau_f; mat(3, 2) = ctau_f; mat(3, 3) = -ctau_f;
    return mat;
}();





// As funções podem usar diretamente essas variáveis

SX S(const casadi::SX& input) {
    // Extrair os elementos do vetor de entrada
    casadi::SX i = input(0);
    casadi::SX j = input(1);
    casadi::SX k = input(2);

    // Criar a matriz de rotação
    casadi::SX out = casadi::SX::zeros(3, 3);
    out(0, 1) = -k;
    out(0, 2) = j;
    out(1, 0) = k;
    out(1, 2) = -i;
    out(2, 0) = -j;
    out(2, 1) = i;

    return out;
}

SX reshape_with_check(const SX& input, size_t new_rows, size_t new_cols, const std::string& name) {
    // Obter tamanho total do input
    size_t total_elements = input.size1() * input.size2();
    size_t new_total_elements = new_rows * new_cols;

    // Adicionar log para depuração
    //std::cout << "[reshape_with_check] Input (" << name << "): rows=" << input.size1()
     //         << ", cols=" << input.size2() << ", total_elements=" << total_elements << std::endl;
    //std::cout << "[reshape_with_check] Novo tamanho solicitado: rows=" << new_rows
        //      << ", cols=" << new_cols << ", total_elements=" << new_total_elements << std::endl;

    // Verificar consistência de dimensões
    if (total_elements != new_total_elements) {
        std::ostringstream oss;
        oss << "Erro em reshape_with_check: Tamanho total inconsistente para reshape de "
            << name << ". Total atual: " << total_elements
            << ", Total solicitado: " << new_total_elements << ".";
        throw std::runtime_error(oss.str());
    }

    // Log para reshape bem-sucedido
    //std::cout << "[reshape_with_check] Sucesso! Realizando reshape." << std::endl;

    // Realizar reshape
    return reshape(input, new_rows, new_cols);
}

SX mtimes_with_check(const SX& A, const SX& B, const std::string& name_A, const std::string& name_B) {
    size_t rows_A = A.size1();
    size_t cols_A = A.size2();
    size_t rows_B = B.size1();
    size_t cols_B = B.size2();

    // Adicionar log para depuração
    //std::cout << "[mtimes_with_check] Matriz A (" << name_A << "): rows=" << rows_A << ", cols=" << cols_A << std::endl;
    //std::cout << "[mtimes_with_check] Matriz B (" << name_B << "): rows=" << rows_B << ", cols=" << cols_B << std::endl;

    // Permitir multiplicação com escalares (1x1)
    if (rows_A == 1 && cols_A == 1) {
        //std::cout << "[mtimes_with_check] Matriz A (" << name_A << ") é um escalar." << std::endl;
        return mtimes(A, B);
    }
    if (rows_B == 1 && cols_B == 1) {
        //std::cout << "[mtimes_with_check] Matriz B (" << name_B << ") é um escalar." << std::endl;
        return mtimes(A, B);
    }

    // Verificar compatibilidade para multiplicação de matrizes
    if (cols_A != rows_B) {
        std::ostringstream oss;
        oss << "Erro em mtimes_with_check: Dimensões incompatíveis para multiplicação entre "
            << name_A << " e " << name_B << ". "
            << "cols(A)=" << cols_A << " deve ser igual a rows(B)=" << rows_B << ".";
        throw std::runtime_error(oss.str());
    }

    // Log para multiplicação bem-sucedida
    //std::cout << "[mtimes_with_check] Sucesso! Realizando multiplicação." << std::endl;

    // Realizar multiplicação
    return mtimes(A, B);
}

SX vertcat_with_check(const std::vector<SX>& components, const std::vector<std::string>& component_names) {
    if (components.empty()) {
        throw std::runtime_error("Nenhum componente fornecido para concatenar.");
    }

    if (components.size() != component_names.size()) {
        throw std::runtime_error("O número de componentes não corresponde ao número de nomes fornecidos.");
    }

    // Verificar dimensões consistentes
    size_t width = components[0].size2();

    // Adicionar log detalhado no início
    //std::cout << "[vertcat_with_check] Verificando " << components.size() << " componentes." << std::endl;

    for (size_t i = 0; i < components.size(); ++i) {
        //std::cout << "[vertcat_with_check] Componente " << i
             //     << " (Nome: " << component_names[i] << "): rows=" << components[i].size1()
                //  << ", cols=" << components[i].size2() << std::endl;

        if (components[i].size2() != width) {
            std::ostringstream oss;
            oss << "Erro em vertcat_with_check: Largura inconsistente detectada no componente " << i
                << " (Nome: " << component_names[i] << "). Esperado: " << width
                << ", Obtido: " << components[i].size2() << ".";
            throw std::runtime_error(oss.str());
        }
    }

    // Exibir uma mensagem antes da concatenação bem-sucedida
    //std::cout << "[vertcat_with_check] Sucesso! Concatenando componentes." << std::endl;

    // Concatenar os componentes verticalmente
    return vertcat(components);
}


SX horzcat_with_check(const std::vector<SX>& components, const std::vector<std::string>& component_names) {
    if (components.empty()) {
        throw std::runtime_error("Nenhum componente fornecido para concatenar.");
    }

    if (components.size() != component_names.size()) {
        throw std::runtime_error("O número de componentes não corresponde ao número de nomes fornecidos.");
    }

    // Verificar consistência nas dimensões verticais (número de linhas)
    size_t height = components[0].size1();

    // Adicionar log detalhado no início
    //std::cout << "[horzcat_with_check] Verificando " << components.size() << " componentes." << std::endl;

    for (size_t i = 0; i < components.size(); ++i) {
        //std::cout << "[horzcat_with_check] Componente " << i
          //        << " (Nome: " << component_names[i] << "): rows=" << components[i].size1()
          //        << ", cols=" << components[i].size2() << std::endl;

        if (components[i].size1() != height) {
            std::ostringstream oss;
            oss << "Erro em horzcat_with_check: Altura inconsistente detectada no componente " << i
                << " (Nome: " << component_names[i] << "). Esperado: " << height
                << ", Obtido: " << components[i].size1() << ".";
            throw std::runtime_error(oss.str());
        }
    }

    // Exibir uma mensagem antes da concatenação bem-sucedida
    //std::cout << "[horzcat_with_check] Sucesso! Concatenando componentes horizontalmente." << std::endl;

    // Concatenar os componentes horizontalmente
    return horzcat(components);
}

//função produto q
SX produto_q(const SX& q1, const SX& q2){
    // Extract components of q1
    SX w = q1(0), x = q1(1), y = q1(2), z = q1(3);
    SX H4_q1 = SX::zeros(4, 4);
    H4_q1(0, 0) = w;  H4_q1(0, 1) = -x; H4_q1(0, 2) = -y; H4_q1(0, 3) = -z;
    H4_q1(1, 0) = x;  H4_q1(1, 1) = w;  H4_q1(1, 2) = -z; H4_q1(1, 3) = y;
    H4_q1(2, 0) = y;  H4_q1(2, 1) = z;  H4_q1(2, 2) = w;  H4_q1(2, 3) = -x;
    H4_q1(3, 0) = z;  H4_q1(3, 1) = -y; H4_q1(3, 2) = x;  H4_q1(3, 3) = w;

    return mtimes_with_check(H4_q1, q2, "H4_q1", "q2");
}

// Funcao verificada
SX mult_qd_casadi(const SX& q1, const SX& q2) {
    // if (q1.size1() != 8 || q1.size2() != 1) {
    //     throw std::invalid_argument("O vetor 'q1' deve ter dimensões 8x1 em mult_qd_casadi.");
    // }
    // if (q2.size1() != 8 || q2.size2() != 1) {
    //     throw std::invalid_argument("O vetor 'q2' deve ter dimensões 8x1 em mult_qd_casadi.");
    // }
    // Separar partes real e dual
    SX r1 = q1(Slice(0, 4)); // Parte real de q1
    SX d1 = q1(Slice(4, 8)); // Parte dual de q1
    SX r2 = q2(Slice(0, 4)); // Parte real de q2
    SX d2 = q2(Slice(4, 8)); // Parte dual de q2

    // Parte real da multiplicação
    SX real_part = SX::vertcat({
        r1(0)*r2(0) - r1(1)*r2(1) - r1(2)*r2(2) - r1(3)*r2(3), // Escalar
        r1(0)*r2(1) + r1(1)*r2(0) + r1(2)*r2(3) - r1(3)*r2(2), // i
        r1(0)*r2(2) - r1(1)*r2(3) + r1(2)*r2(0) + r1(3)*r2(1), // j
        r1(0)*r2(3) + r1(1)*r2(2) - r1(2)*r2(1) + r1(3)*r2(0)  // k
    });

    // Parte dual da multiplicação
    SX dual_part = SX::vertcat({
        // Escalar
        r1(0)*d2(0) + d1(0)*r2(0) - r1(1)*d2(1) - r1(2)*d2(2) - r1(3)*d2(3)
        - d1(1)*r2(1) - d1(2)*r2(2) - d1(3)*r2(3),
        // i
        r1(0)*d2(1) + r1(1)*d2(0) + r1(2)*d2(3) - r1(3)*d2(2)
        + d1(0)*r2(1) + d1(1)*r2(0) + d1(2)*r2(3) - d1(3)*r2(2),
        // j
        r1(0)*d2(2) - r1(1)*d2(3) + r1(2)*d2(0) + r1(3)*d2(1)
        + d1(0)*r2(2) - d1(1)*r2(3) + d1(2)*r2(0) + d1(3)*r2(1),
        // k
        r1(0)*d2(3) + r1(1)*d2(2) - r1(2)*d2(1) + r1(3)*d2(0)
        + d1(0)*r2(3) + d1(1)*r2(2) - d1(2)*r2(1) + d1(3)*r2(0)
    });

    // Concatenar parte real e dual
    SX result = SX::vertcat({real_part, dual_part});
    return result;
}

SX primario_qd_casadi(const SX& q) {
    SX q_aux = SX::vertcat({q(0), q(1), q(2), q(3), SX(0), SX(0), SX(0), SX(0)});
    return q_aux;
}

SX dual_qd_casadi(const SX& q) {
    SX q_aux = SX::vertcat({q(4), q(5), q(6), q(7), SX(0), SX(0), SX(0), SX(0)});
    return q_aux;
}

// Funcao verificada
SX conj_qd_casadi(const casadi::SX& q){
    // Verifique se q é um vetor linha com 8 elementos
    if (q.size1() != 8 || q.size2() != 1) {
        std::cout << "q tem tamanho errado: " << q << std::endl;
        throw std::invalid_argument("O quaternio dual 'q' deve ser um vetor linha com 8 componentes (8x1) conjugado.");
    }

    // Componentes individuais do quaternio dual
    SX q0 = q(0);      // Componente escalar
    SX q1 = -q(1);     // Componente x negada
    SX q2 = -q(2);     // Componente y negada
    SX q3 = -q(3);     // Componente z negada
    SX q4 = q(4);      // Componente dual escalar
    SX q5 = -q(5);     // Componente x dual negada
    SX q6 = -q(6);     // Componente y dual negada
    SX q7 = -q(7);     // Componente z dual negada

    // Formar o quaternio conjugado como vetor linha
    SX q_conj = SX::vertcat({q0, q1, q2, q3, q4, q5, q6, q7});

    // Garantir que a saída também seja um vetor linha
    if (q_conj.size1() != 8 || q_conj.size2() != 1) {
        throw std::runtime_error("Erro interno: q_inv não é um vetor linha (8x1) conjugado.");
    }
    return q_conj;
}

// Funcao verificada
SX ad_qd_casadi(const SX& q1, const SX& q2)
    {
    // Calcular o conjugado de q1
    SX q1_conj = conj_qd_casadi(q1);

    // Primeira Parte --> q1* q2
    SX p1 = mult_qd_casadi(q1, q2);

    // Segunda Parte --> p1 * q1_conjugado
    SX ad_qd = mult_qd_casadi(p1, q1_conj);

    return ad_qd;
}

// Funcao verificada
// Extrai a parte vetorial (x, y, z) da parte real de um quaternion dual
SX vec3_qd_casadi(const SX& q) {
    // Assume q é um vetor linha ou coluna com 8 elementos
    // Parte real: q(0), q(1), q(2), q(3)
    // Parte vetorial: q(1), q(2), q(3)
    return SX::vertcat({q(1), q(2), q(3)});
}

// Funcao verificada
// Implementa M_qd_casadi equivalente à função MATLAB fornecida
SX m_qd_casadi(const SX& m, const SX& q) {
    SX q1 = m * vec3_qd_casadi(q); // Multiplica um escalar ou uma matriz pela parte vetorial
    // Monta o vetor resultado conforme especificado
    return SX::vertcat({SX(0), q1(0), q1(1), q1(2), SX(0), SX(0), SX(0), SX(0)});
}

// Funcao verificada
// Soma de quaterniões duais compatível com CasADi
SX soma_qd_casadi(const SX& q1, const SX& q2) {
    // Soma das partes
    SX real_vec = q1 + q2;

    return real_vec;
}

// Funcao verificada
// Subtração de quaterniões duais compatível com CasADi
SX sub_qd_casadi(const SX& q1, const SX& q2) {
    // Subtração das partes
    SX real_vec = q1 - q2;

    return real_vec;
}

// Funcao verificada
// Produto vetorial de quaternios duais compatível com CasADi
SX cross_qd_casadi(const SX& q1, const SX& q2) {
    // Calcula q1*q2 e q2*q1
    SX parte_1 = mult_qd_casadi(q1, q2);
    SX parte_2 = mult_qd_casadi(q2, q1);
    // Subtração
    SX sub = parte_1 - parte_2;
    // Separar parte real e parte dual
    SX real_part = sub(Slice(0, 4));
    SX dual_part = sub(Slice(4, 8));
    // Multiplicar por 0.5
    real_part = real_part * 0.5;
    dual_part = dual_part * 0.5;
    // Retornar combinação das partes real e dual
    return SX::vertcat({real_part, dual_part});
}

// Funcao verificada
SX p_qd_casadi(const SX& q) {
    // Extrai parte real (r) e dual (d)
    SX r = q(Slice(0, 4));   // Parte real
    SX d = q(Slice(4, 8));   // Parte dual

    // Inverso do quaternion r (apenas conjugado, pois r é unitário)
    SX r_inv = SX::vertcat({r(0), -r(1), -r(2), -r(3)});

    // produto_q espera vetores linha (1x4), então transpor d e r_inv se necessário
    SX p_quat = 2 * produto_q(d, r_inv);

    // A parte vetorial de p_quat é o vetor de posição p (índices 1,2,3)
    SX p = SX::horzcat({p_quat(1), p_quat(2), p_quat(3)});
    return p;
}

// Funcao verificada
// Norma de quaternio dual compatível com CasADi
SX norma_qd_casadi(const SX& q) {
    // Quaternio dual conjugado
    SX q_conj = conj_qd_casadi(q);

    // Multiplicação de q pelo seu conjugado
    SX m = mult_qd_casadi(q, q_conj);

    SX r1 = m(Slice(0, 4));
    SX d1 = m(Slice(4, 8));

    // Norma real (escalar)
    SX r1_scalar = sqrt(r1(0));
    // Norma dual (escalar)
    SX d1_scalar = d1(0) / (2 * r1_scalar);

    // Resultado: vetor linha 1x2 [norma_real, norma_dual]
    SX res = SX::horzcat({r1_scalar, d1_scalar});
    return res;
}

SX mtimes_qd_casadi(const SX& m, const SX& q) {
    // Chama a função vec3_qd_casadi (assumindo que existe)
    if (m.size1() != 3 || m.size2() != 3) {
        throw std::invalid_argument("A matriz 'm' deve ter dimensões 3x3 em mtimes_qd_casadi.");
    }
    if (q.size1() != 8 || q.size2() != 1) {
        throw std::invalid_argument("O vetor 'q' deve ter dimensões 8x1 em mtimes_qd_casadi.");
    }

    SX q1 = m(0,0) * q(1);
    // std::cout << "q1, m, q: " << q1 << ' ' <<m(0,0) << ' ' << q(1) << std::endl; //pode excluir
    SX q2 = m(1,1) * q(2);
    // std::cout << "q2: " << q2 << std::endl; //pode excluir
    SX q3 = m(2,2) * q(3);
    // std::cout << "q3: " << q3 << std::endl; //pode excluir
    
    // Constrói o vetor resultado
    SX res = SX::vertcat({
        SX(0),           // elemento 1: 0
        q1,           // elemento 2: q1(1) - note: índices começam em 0
        q2,           // elemento 3: q1(2)
        q3,           // elemento 4: q1(3)
        SX(0), SX(0), SX(0), SX(0)  // elementos 5-8: zeros
    });
    
    return res;
}

// // Funcao problematica
// SX normalize_qd_casadi(const SX& qd) {
//     if (qd.size1() != 8 || qd.size2() != 1) {
//         throw std::invalid_argument("O vetor 'qd' deve ter dimensões 8x1 em normalize.");
//     }
//     // Separar partes real e dual
//     SX q_real = qd(Slice(0, 4));
//     SX q_dual = qd(Slice(4, 8));

//     // Calcular norma da parte real
//     SX norm_real = sqrt(mtimes(q_real, q_real.T()));

//     // Caso especial: norma muito pequena
//     SX small_norm = norm_real < 1e-12;
//     SX scale_real = if_else(small_norm, SX(1), SX(1) / norm_real);
//     SX q_real_norm = scale_real * q_real;

//     // Corrigir parte dual (garantir perpendicularidade)
//     SX dot_prod = mtimes(q_real_norm, q_dual);
//     SX q_dual_norm = (q_dual - dot_prod * q_real_norm) * scale_real;

//     // Definir quatérnio identidade dual caso necessário
//     SX q_real_norm_safe = if_else(norm_real < small_norm, SX::vertcat({1, 0, 0, 0}), q_real_norm);
//     SX q_dual_norm_safe = if_else(norm_real < small_norm, SX::vertcat({0, 0, 0, 0}), q_dual_norm);

//     // Concatenar resultado (vetor linha 1x8)
//     SX qd_norm = SX::vertcat({q_real_norm_safe, q_dual_norm_safe});
//     return qd_norm;
// }

// Funcao verificada
// Exponencial de quaternio dual compatível com CasADi
SX exp_qd_casadi(const SX& q) {
    SX n = norma_qd_casadi(q); // n(0): norma real, n(1): norma dual
    SX small_value = 1e-10;
    SX n1 = n(0) + small_value * (n(0) == 0);

    // p_exp = cos(n(0)) * [1;0;0;0] + (sin(n(0))/n1) * q(0:4)
    SX p_exp = cos(n(0)) * SX::vertcat({1, 0, 0, 0}) + (sin(n(0)) / n1) * q(Slice(0, 4));

    // Criação dos vetores real e dual (8 elementos cada)
    SX real = SX::vertcat({p_exp, SX::zeros(4, 1)});
    // std::cout << "real: " << real << std::endl; //pode excluir
    SX dual = SX::vertcat({q(Slice(4, 8)), SX::zeros(4, 1)});
    // std::cout << "dual: " << dual << std::endl; //pode excluir

    // Multiplicação do dual pelo real
    SX ex_dual = mult_qd_casadi(dual, real);

    // Resultado final: [ex_real; ex_dual(0:3)]
    SX res = SX::vertcat({p_exp, ex_dual(Slice(0, 4))});
    return res;
}


// Funcao verificada
// Inversa de quaternio dual compatível com CasADi (entrada e saída: SX linha 1x8)
SX inv_qd_casadi(const SX& q) {
    // Quaternio dual conjugado
    SX q_conj = conj_qd_casadi(q);

    // Multiplicação de q pelo seu conjugado (norma dual ao quadrado)
    SX x = mult_qd_casadi(q, q_conj);

    // Extração das partes real e dual do resultado
    SX xr = x(Slice(0, 4));  // Parte real
    SX xd = x(Slice(4, 8));  // Parte dual

    // Parte escalar (norma real e dual)
    SX r = xr(0);  // Parte real escalar
    SX d = xd(0);  // Parte dual escalar

    // Evitar divisão por zero
    SX small_value = 1e-10;
    SX r_safe = r + small_value * (r == 0);

    // Cálculo da inversa (expansão de Taylor)
    SX inv_primary = 1 / r_safe;
    SX inv_dual = -d / (r_safe * r_safe);

    // Cálculo das partes real e dual da inversa
    SX real_part_inv = q_conj(Slice(0, 4)) * inv_primary;
    SX dual_part_inv = q_conj(Slice(4, 8)) * inv_primary + q_conj(Slice(0, 4)) * inv_dual;

    // Resultado final: vetor linha 1x8
    SX res = SX::vertcat({real_part_inv, dual_part_inv});
    return res;
}


// Logaritmo de quaternio dual compatível com CasADi
SX log_qd_casadi(const SX& dq) {
    if (dq.size1() != 8 || dq.size2() != 1) {
        throw std::invalid_argument("O vetor 'dq' deve ter dimensões 8x1 em log.");
    }
    // Extrair parte real (r) e dual (d)
    SX r = dq(Slice(0, 4));
    SX d = dq(Slice(4, 8));

    // Normalizar a parte real
    SX norm_r = sqrt(mtimes(r.T(), r));
    SX r_norm = r / norm_r;

    // Componentes da parte real
    SX s = r_norm(0);
    SX v = SX::vertcat({r_norm(1), r_norm(2), r_norm(3)});
    SX nv_sq = mtimes(v.T(), v); // ||v||
    SX nv = sqrt(nv_sq);         // ||v||^2
    // std::cout << "s: " << s << ", v: " << v << ", nv_sq: " << nv_sq << ", nv: " << nv << std::endl; //pode excluir

    // Cálculo seguro de theta (evita divisão por zero)
    SX epsilon = 1e-12;
    SX safe_nv = if_else(nv > epsilon, nv, epsilon);
    SX theta = 2 * atan2(safe_nv, s);
    // std::cout << "theta: " << theta << std::endl; //pode excluir

    // Cálculo do logaritmo real (sem singularidade)
    SX scale = if_else(nv_sq > epsilon, theta / safe_nv, 0);
    SX log_real = scale * v;
    // std::cout << "log_real: " << log_real << std::endl; //pode excluir

    // Cálculo do logaritmo dual
    SX r_conj = SX::vertcat({r_norm(0), -r_norm(1), -r_norm(2), -r_norm(3)});
    SX q_temp = produto_q(r_conj, d);
    // std::cout << "log_dual before real: " << log_real<< std::endl; //pode excluir
    SX log_dual = SX::vertcat({q_temp(1), q_temp(2), q_temp(3)}); // Parte vetorial
    // std::cout << "log_dual before scale: " << r_conj << std::endl; //pode excluir


    SX res = SX::vertcat({SX(0), log_real, SX(0), log_dual});
    // std::cout << "log_dual before dps: " << log_dual << std::endl;  //pode excluir
    return res;
}

// Funcao verificada
SX rot_qd_casadi(const SX& q) {
    if (q.size1() != 8 || q.size2() != 1) {
        throw std::invalid_argument("O vetor 'q' deve ter dimensões 8x1 em rot_qd.");
    }
    // Parte real: q(0:3), parte dual: zeros
    SX result = SX::vertcat({q(Slice(0, 4)), SX::zeros(4, 1)});
    // std::cout << "rot_qd_casadi result: " << result << std::endl; //pode excluir
    return result;
}

// Funcao verificada
SX translacao_qd_casadi(const SX& q) {
    if (q.size1() != 8 || q.size2() != 1) {
        throw std::invalid_argument("O vetor 'q' deve ter dimensões 8x1 em translacao.");
    }
    // Separar partes real e dual
    SX real_part = q(Slice(0, 4));  // Parte real
    SX dual_part = q(Slice(4, 8));  // Parte dual

    // Conjugado da parte real
    SX real_conj = SX::vertcat({real_part(0), -real_part(1), -real_part(2), -real_part(3)});

    // Cria quaternios QD para a parte real e dual (1x8)
    SX q_real = SX::vertcat({real_conj, SX::zeros(4, 1)});
    SX q_dual = SX::vertcat({dual_part, SX::zeros(4, 1)});

    // Multiplica a parte dual pela parte real conjugada e escala por 2
    SX p = mult_qd_casadi(q_dual, q_real);
    p = soma_qd_casadi(p, p); // Escalar por 2

    // Retorna o resultado como quaternio dual (1x8)
    return p;
}


SX verify_and_normalize_quaternion(const SX& q, double tol = 1e-6) {
    // Calcula a norma do quaternion
    SX norm_q = norm_2(q);

    // Verifica se a norma é zero (evita divisão por zero)
    if (DM(norm_q).is_zero()) {
        std::cerr << "Erro: Quaternion inválido com norma zero." << std::endl;
        throw std::runtime_error("Norma do quaternion é zero.");
    }

    // Verifica se o quaternion é unitário dentro da tolerância
    SX diff = fabs(norm_q - SX(1));
    if (DM(diff).scalar() > tol) {
        std::cerr << "Erro: Quaternion não é unitário. Norma = " << DM(norm_q) << std::endl;
        throw std::runtime_error("Norma do quaternion não está dentro da tolerância para ser considerado unitário.");
    }

    // Retorna o quaternion normalizado
    return q / norm_q;
}

//transformar euler pra quaternion
Eigen::Quaterniond eul2quat(double roll, double pitch, double yaw)
{
// Create individual rotations using Eigen's AngleAxis
Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());  // Rotation around X-axis
Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY()); // Rotation around Y-axis
Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());    // Rotation around Z-axis

// Combine the rotations in the order: yaw, pitch, roll (common aerospace convention)
Eigen::Quaterniond q = yawAngle * pitchAngle * rollAngle;

// Return the resulting quaternion
return q;
}


SX quaternion_inverse(const casadi::SX& q) {
    // Verifique se q é um vetor linha com 4 elementos
    if (q.size1() != 1 || q.size2() != 4) {
        throw std::invalid_argument("O quaternio 'q' deve ser um vetor linha com 4 componentes (1x4).");
    }

    // Norma do quaternio
    SX norm_q = norm_2(q);

    // Componentes individuais do quaternio
    SX q0 = q(0);      // Componente escalar
    SX q1 = -q(1);     // Componente x negada
    SX q2 = -q(2);     // Componente y negada
    SX q3 = -q(3);     // Componente z negada

    // Formar o quaternio inverso como vetor linha
    SX q_inv = SX::horzcat({q0, q1, q2, q3}) / (norm_q * norm_q);

    // Garantir que a saída também seja um vetor linha
    if (q_inv.size1() != 1 || q_inv.size2() != 4) {
        throw std::runtime_error("Erro interno: q_inv não é um vetor linha (1x4).");
    }

    return q_inv;
}


// Function to multiply two quaternions
SX quatmultiply(const SX& q1, const SX& q2) {
    // Verificar se q1 e q2 são vetores linha com 4 componentes
    if (q1.size1() != 1 || q1.size2() != 4) {
        throw std::invalid_argument("O quaternio 'q1' deve ser um vetor linha com 4 componentes (1x4).");
    }
    if (q2.size1() != 1 || q2.size2() != 4) {
        throw std::invalid_argument("O quaternio 'q2' deve ser um vetor linha com 4 componentes (1x4).");
    }

    // Extrair os componentes de q1
    SX w1 = q1(0), x1 = q1(1), y1 = q1(2), z1 = q1(3);

    // Extrair os componentes de q2
    SX w2 = q2(0), x2 = q2(1), y2 = q2(2), z2 = q2(3);

    // Fórmula de multiplicação de quaterniões
    SX w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2;
    SX x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2;
    SX y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2;
    SX z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2;

    // Combinar os componentes em um único quaternio (vetor linha)
    SX q_result_mult = SX::horzcat({w, x, y, z});

    // Garantir que a saída também seja um vetor linha (1x4)
    if (q_result_mult.size1() != 1 || q_result_mult.size2() != 4) {
        throw std::runtime_error("Erro interno: o quaternio resultante não é um vetor linha (1x4).");
    }

    return q_result_mult;
}



// Function to calculate the orientation error between two quaternions
SX quaternion_error(const SX& q, const SX& qd) {
    // Verificar se q e qd são vetores linha com 4 componentes
    if (q.size1() != 1 || q.size2() != 4) {
        throw std::invalid_argument("O quaternio 'q' deve ser um vetor linha com 4 componentes (1x4).");
    }
    if (qd.size1() != 1 || qd.size2() != 4) {
        throw std::invalid_argument("O quaternio 'qd' deve ser um vetor linha com 4 componentes (1x4).");
    }

    // qtil = q * qd^(-1)
    SX q_inv = quaternion_inverse(q);  // Inverse of q
    SX qtil = quatmultiply(q_inv, qd); // Multiply q_inv with qd

    // Garantir que qtil seja um vetor linha (1x4)
    if (qtil.size1() != 1 || qtil.size2() != 4) {
        throw std::runtime_error("Erro interno: o quaternio de erro 'qtil' não é um vetor linha (1x4).");
    }

    return qtil;
}


SX compute_Rq(const SX& X_k_aux) {
    // Inicializa a matriz Rq com dimensões 3x3 preenchida com zeros
    SX Rq = SX::zeros(3, 3);

    // Preenche os elementos da matriz de rotação Rq
    Rq(0, 0) = 1 - 2 * pow(X_k_aux(8), 2) - 2 * pow(X_k_aux(9), 2);
    Rq(0, 1) = 2 * X_k_aux(7) * X_k_aux(8) - 2 * X_k_aux(9) * X_k_aux(6);
    Rq(0, 2) = 2 * X_k_aux(7) * X_k_aux(9) + 2 * X_k_aux(8) * X_k_aux(6);

    Rq(1, 0) = 2 * X_k_aux(7) * X_k_aux(8) + 2 * X_k_aux(9) * X_k_aux(6);
    Rq(1, 1) = 1 - 2 * pow(X_k_aux(7), 2) - 2 * pow(X_k_aux(9), 2);
    Rq(1, 2) = 2 * X_k_aux(8) * X_k_aux(9) - 2 * X_k_aux(7) * X_k_aux(6);

    Rq(2, 0) = 2 * X_k_aux(7) * X_k_aux(9) - 2 * X_k_aux(8) * X_k_aux(6);
    Rq(2, 1) = 2 * X_k_aux(8) * X_k_aux(9) + 2 * X_k_aux(7) * X_k_aux(6);
    Rq(2, 2) = 1 - 2 * pow(X_k_aux(7), 2) - 2 * pow(X_k_aux(8), 2);

    return Rq;
}

SX compute_cross1(const SX& ohm_a, const SX& j_ohm) {
    // Verificar se os vetores têm dimensão 3x1
    if (ohm_a.size1() != 3 || ohm_a.size2() != 1) {
        throw std::invalid_argument("O vetor 'ohm_a' deve ter dimensões 3x1.");
    }
    if (j_ohm.size1() != 3 || j_ohm.size2() != 1) {
        throw std::invalid_argument("O vetor 'j_ohm' deve ter dimensões 3x1.");
    }

    // Inverter os componentes de ohm_a (multiplicação por -1)
    SX neg_ohm_a = -ohm_a;

    // Produto vetorial entre neg_ohm_a e j_ohm
    SX cross1 = SX::vertcat({
        neg_ohm_a(1) * j_ohm(2) - neg_ohm_a(2) * j_ohm(1), // Componente x
        neg_ohm_a(2) * j_ohm(0) - neg_ohm_a(0) * j_ohm(2), // Componente y
        neg_ohm_a(0) * j_ohm(1) - neg_ohm_a(1) * j_ohm(0)  // Componente z
    });

    return cross1;
}

SX compute_psiq(const SX& q_a, const SX& P_extra) {
    // Obtém a orientação desejada a partir de P_extra
    SX qd = P_extra(Slice(19, 23)); // Obtém a orientação desejada (sem orientação com relação ao frame inercial)

    SX aux = SX::horzcat({1, 0, 0, 0});

    // Calcula o erro de quaternion
    SX qtil = aux - quaternion_error(q_a.T(), qd.T());

    // Calcula psiq como o produto interno
    SX psiq = mtimes_with_check(qtil, qtil.T(), "qtil", "qtil.T()");

    return psiq;
}

SX compute_xbarra_qd(const SX& pose, const SX& P) {
    // Calcula os elementos de xbarra
    SX unit = SX::vertcat({SX(1), SX(0), SX(0), SX(0), SX(0), SX(0), SX(0), SX(0)});
    SX q_e_hat = mult_qd_casadi(inv_qd_casadi(P(Slice(16, 24), 0)), pose);
     SX x_till = unit - q_e_hat;
    //SX x_till = log_qd_casadi(q_e_hat);

    return x_till;
}

SX compute_hbarra_qd(const SX& pose, const SX& heligiro, const SX& P) {
    // Calcula os elementos de hbarra
    SX omega_b = heligiro(Slice(1, 4)); // Vetor ω_b (3x1)
    SX p_dot_b = heligiro(Slice(5, 8)); // Vetor ṗ_b (3x1)
    SX p_b = pose(Slice(5, 8));       // Vetor p_b (3x1)

    SX termo_coriolis = compute_cross1(omega_b, p_b); // ω_b × p_b (3x1)
    SX parte_dual_xi_b = p_dot_b + termo_coriolis;  // ṗ_b + ω_b × p_b (3x1)

    SX heligiro_primario = heligiro(Slice(0, 4));
    SX xi_b = SX::vertcat({heligiro_primario,  // Parte real: ω_b (4x1 - quaternion)
                         SX(0), 
                         parte_dual_xi_b});        // Vetor da parte dual (3x1)
    
    SX q_e_hat = mult_qd_casadi(inv_qd_casadi(P(Slice(16, 24), 0)), pose);
    SX xi_d = ad_qd_casadi(conj_qd_casadi(q_e_hat), P(Slice(24, 32), 0));
    SX erro_heligiro = xi_b - xi_d;

    return erro_heligiro;
}

SX compute_hbarra_qd2(const SX& pose, const SX& heligiro, const SX& P) {
    // Calcula os elementos de hbarra
    SX omega_i = SX::vertcat({SX(0), heligiro(Slice(1, 4)), SX::zeros(4,1)}); // Vetor ω_i (4x1 - quaternion)
    SX omega_b = ad_qd_casadi(conj_qd_casadi(rot_qd_casadi(pose)), omega_i); // Vetor ω_b (3x1)
    omega_b = omega_b(Slice(1, 4)); // Extrai apenas a velocidade angular (3x1)

    SX p_i = translacao_qd_casadi(pose); // Vetor p_i (8x1)
    SX p_b = ad_qd_casadi(conj_qd_casadi(rot_qd_casadi(pose)), p_i);       // Vetor p_b (3x1)
    p_b = p_b(Slice(1, 4));     // Extrai apenas a parte vetorial (3x1)

    SX p_dot_b_auxiliar = ad_qd_casadi(conj_qd_casadi(pose), heligiro); // Vetor ṗ_i (8x1)
    SX p_dot_b = p_dot_b_auxiliar(Slice(5,8));         // Vetor ṗ_b (3x1)

    SX termo_coriolis = compute_cross1(omega_b, p_b); // ω_b × p_b (3x1)
    SX parte_dual_xi_b = p_dot_b + termo_coriolis;  // ṗ_b + ω_b × p_b (3x1)

    SX heligiro_primario = heligiro(Slice(0, 4));
    SX xi_b = SX::vertcat({heligiro_primario,  // Parte real: ω_b (4x1 - quaternion)
                         SX(0), 
                         parte_dual_xi_b});        // Vetor da parte dual (3x1)
    
    SX q_e_hat = mult_qd_casadi(inv_qd_casadi(P(Slice(16, 24), 0)), pose);
    SX xi_d = ad_qd_casadi(conj_qd_casadi(q_e_hat), P(Slice(24, 32), 0));
    SX erro_heligiro = xi_b - xi_d;

    return erro_heligiro;
}


SX compute_J_qd(const SX& Q, const SX& Qh, const SX& R, const SX& X0, const SX& X_r, const SX& X_k, const SX& U_k) {
    // Cálculo de P
    SX P = vertcat_with_check({X0, X_r}, {"X0", "X_r"});

    // Cálculo de xbarra
    SX xbarra = compute_xbarra_qd(X_k(Slice(0, 8)), P);
    SX hbarra = compute_hbarra_qd2(X_k(Slice(0, 8)), X_k(Slice(8, 16)), P);

    SX conJ = mtimes(mat_aux_con, U_k); //con
    SX torque_motorJ = conJ(0);
    SX tauJ = SX::vertcat({conJ(1), conJ(2), conJ(3)});

    // Cálculo de ubarra
    SX ubarra = SX::zeros(4, 1);
    ubarra(0, 0) = torque_motorJ - m * gr;         // Primeiro componente fi
    ubarra(1, 0) = tauJ(0);          // Primeiro componente de tau
    ubarra(2, 0) = tauJ(1);          // Segundo componente de tau
    ubarra(3, 0) = tauJ(2);          // Terceiro componente de tau    

    // Cálculo incremental do custo J
    SX J =    mtimes_with_check(xbarra.T(), mtimes_with_check(Q, xbarra, "Q", "xbarra"), "xbarra.T()", "Q * xbarra")
            + mtimes_with_check(hbarra.T(), mtimes_with_check(Qh, hbarra, "Qh", "hbarra"), "hbarra.T()", "Qh * hbarra")
            + mtimes_with_check(ubarra.T(), mtimes_with_check(R, ubarra, "R", "ubarra"), "ubarra.T()", "R * ubarra");

    return J;
}

SX compute_J_qd_rapido(const SX& Q, const SX& Qh, const SX& R, const SX& X0, const SX& X_r, const SX& X_k, const SX& U_k, int k, const float zeta) {
    // Cálculo de P
    SX P = vertcat_with_check({X0, X_r}, {"X0", "X_r"});

    // Cálculo de xbarra
    SX xbarra = compute_xbarra_qd(X_k(Slice(0, 8)), P);
    SX hbarra = compute_hbarra_qd2(X_k(Slice(0, 8)), X_k(Slice(8, 16)), P);

    SX conJ = mtimes(mat_aux_con, U_k); //con
    SX torque_motorJ = conJ(0);
    SX tauJ = SX::vertcat({conJ(1), conJ(2), conJ(3)});

    // Cálculo de ubarra
    SX ubarra = SX::zeros(4, 1);
    ubarra(0, 0) = torque_motorJ - m * gr;         // Primeiro componente fi
    ubarra(1, 0) = tauJ(0);          // Primeiro componente de tau
    ubarra(2, 0) = tauJ(1);          // Segundo componente de tau
    ubarra(3, 0) = tauJ(2);          // Terceiro componente de tau

    // Cálculo incremental do custo J
    SX J =    (mtimes_with_check(xbarra.T(), mtimes_with_check(Q, xbarra, "Q", "xbarra"), "xbarra.T()", "Q * xbarra"))*std::pow(zeta, k)
            + mtimes_with_check(hbarra.T(), mtimes_with_check(Qh, hbarra, "Qh", "hbarra"), "hbarra.T()", "Qh * hbarra")
            + mtimes_with_check(ubarra.T(), mtimes_with_check(R, ubarra, "R", "ubarra"), "ubarra.T()", "R * ubarra");

    return J;
}

std::tuple<SX, SX> compute_next_qd(const SX& state, const SX& U, SX& aak, const double amost) {
    SX gravidade = SX::vertcat({SX(0), SX(0), SX(0), -9.81}); // Parte real com zero na parte dual

    // Criar a matriz de inércia
    SX J_inertia_sx = SX::zeros(3, 3);
    J_inertia_sx(0,0) = Ixx;
    J_inertia_sx(1,1) = Iyy;
    J_inertia_sx(2,2) = Izz;

    SX J_inertia_inv_sx = SX::zeros(3, 3);
    J_inertia_inv_sx(0,0) = Ixxinv;
    J_inertia_inv_sx(1,1) = Iyyinv;
    J_inertia_inv_sx(2,2) = Izzinv;

    SX pose_i = state(Slice(0, 8));
    SX heligiro_i = state(Slice(8, 16));

    SX w_i = ad_qd_casadi(rot_qd_casadi(conj_qd_casadi(pose_i)), primario_qd_casadi(heligiro_i));
    SX v_i = dual_qd_casadi(heligiro_i) - cross_qd_casadi(translacao_qd_casadi(pose_i), primario_qd_casadi(heligiro_i));

    SX con = mtimes(mat_aux_con, U); //con
    SX tauk = SX::vertcat({SX(0), con(1), con(2), con(3), SX::zeros(4,1)}); //torque
    SX ftk = SX::horzcat({SX(0), SX(0), SX(0), con(0), SX::zeros(1,4)}); //empuxo

    pose_i = mult_qd_casadi(exp_qd_casadi(amost*0.5*heligiro_i), pose_i);
    // std::cout << "pose_i depois: " << pose_i << std::endl;
    SX af = (1/m)*ad_qd_casadi(rot_qd_casadi(pose_i), ftk);
    // std::cout << "af: " << af << std::endl;
    SX at = cross_qd_casadi(v_i, primario_qd_casadi(heligiro_i)) + cross_qd_casadi(translacao_qd_casadi(pose_i), ad_qd_casadi(rot_qd_casadi(pose_i), aak));
    // std::cout << "at: " << at << std::endl;
    SX aak_auxiliar = ad_qd_casadi(rot_qd_casadi(pose_i), aak);
    heligiro_i = heligiro_i + amost* SX::vertcat({aak_auxiliar(Slice(0,4)), (af(Slice(0,4)) + at(Slice(0,4)) + gravidade)});
    SX M1 = mtimes_qd_casadi(J_inertia_sx, w_i);

    SX cross_p = cross_qd_casadi(-w_i, M1);

    SX soma = cross_p + tauk;

    aak = mtimes_qd_casadi(J_inertia_inv_sx, soma);

    SX estados_t = SX::vertcat({pose_i, heligiro_i});
    // std::cout << "estados_t: " << estados_t << std::endl;

    return std::make_tuple(estados_t, aak);
}

SX compute_xbarra(const SX& e_a, const SX& v_a, const SX& ohm_a, const SX& P_extra, const SX& psiq) {
    // Inicializa o vetor xbarra com 10 elementos, todos zerados
    SX xbarra = SX::zeros(10, 1);

    // Calcula os elementos de xbarra usando as diferenças e valores fornecidos
    xbarra(0, 0) = e_a(0, 0) - P_extra(13, 0);  // Primeiro elemento
    xbarra(1, 0) = e_a(1, 0) - P_extra(14, 0);  // Segundo elemento
    xbarra(2, 0) = e_a(2, 0) - P_extra(15, 0);  // Terceiro elemento
    xbarra(3, 0) = v_a(0, 0) - P_extra(16, 0);  // Quarto elemento
    xbarra(4, 0) = v_a(1, 0) - P_extra(17, 0);  // Quinto elemento
    xbarra(5, 0) = v_a(2, 0) - P_extra(18, 0);  // Sexto elemento
    xbarra(6, 0) = psiq;                        // Sétimo elemento
    xbarra(7, 0) = ohm_a(0, 0) - P_extra(23, 0); // Oitavo elemento
    xbarra(8, 0) = ohm_a(1, 0) - P_extra(24, 0); // Nono elemento
    xbarra(9, 0) = ohm_a(2, 0) - P_extra(25, 0); // Décimo elemento

    return xbarra;
}

SX compute_J(const SX& Q, const SX& R, const SX& X0, const SX& X_r, const SX& X_k, const SX& U_k) {
    // Cálculo de P_extra
    SX P_extra = vertcat_with_check({X0, X_r}, {"X0", "X_r"});

    // Cálculo de psiq
    SX psiq = compute_psiq(X_k(Slice(6, 10)), P_extra);

    // Cálculo de xbarra
    SX xbarra = compute_xbarra(X_k(Slice(0, 3)), X_k(Slice(3, 6)), X_k(Slice(10, 13)), P_extra, psiq);

    SX conJ = mtimes(mat_aux_con, U_k); //con
    SX torque_motorJ = conJ(0);
    //SX T = SX::vertcat({SX(0), SX(0), conJ(0)}); //T
    SX tauJ = SX::vertcat({conJ(1), conJ(2), conJ(3)});

    // Cálculo da força resultante fi
    SX fi = mtimes_with_check(compute_Rq(X_k.T()), e_z, "Rq", "e_z") * torque_motorJ - m * gr * e_z;



    // Cálculo de ubarra
    SX ubarra = SX::zeros(4, 1);
    ubarra(0, 0) = norm_2(fi);         // Primeiro componente fi
    ubarra(1, 0) = tauJ(0);          // Primeiro componente de tau
    ubarra(2, 0) = tauJ(1);          // Segundo componente de tau
    ubarra(3, 0) = tauJ(2);          // Terceiro componente de tau

    // Cálculo incremental do custo J
    SX J = mtimes_with_check(xbarra.T(), mtimes_with_check(Q, xbarra, "Q", "xbarra"), "xbarra.T()", "Q * xbarra")
            + mtimes_with_check(ubarra.T(), mtimes_with_check(R, ubarra, "R", "ubarra"), "ubarra.T()", "R * ubarra");// + mtimes(compute_psiq(X_k(Slice(6, 10)), P_extra),Kpsi;, 200, "psiq", "200");

    return J;
}

SX compute_next_euler(const SX& state, const SX& u1, const double amost) {

    // Criar a matriz de inércia
    SX J_inertia_sx = SX::zeros(3, 3);
    J_inertia_sx(0,0) = Ixx;
    J_inertia_sx(1,1) = Iyy;
    J_inertia_sx(2,2) = Izz;

    SX J_inertia_inv_sx = SX::zeros(3, 3);
    J_inertia_inv_sx(0,0) = Ixxinv;
    J_inertia_inv_sx(1,1) = Iyyinv;
    J_inertia_inv_sx(2,2) = Izzinv;

    // Estados atuais
    SX e_a = state(Slice(0, 3));  // Posição
    SX v_a = state(Slice(3, 6));    // Velocidades
    SX q_a = state(Slice(6,10));
    SX ohm_a = state(Slice(10, 13)); // Velocidades angulares


    // Calculando as forças e torques
    SX con = mtimes(mat_aux_con, u1); //con
    SX torque_motor = con(0);
    SX T = SX::vertcat({SX(0), SX(0), con(0)}); //T
    SX tau = SX::vertcat({con(1), con(2), con(3)}); //tau

    // Matriz de rotação (Rq)
    SX Rq = compute_Rq(state.T());

    // Acelerações lineares
    SX Ident = SX::eye(3);
    SX j_ohm = mtimes_with_check(J_inertia_sx, ohm_a, "J_inertia_sx", "ohm_a");
    //SX cross1 = mtimes_with_check(S(-ohm_a),j_ohm, "S(-ohm_a)", "j_ohm");
    SX cross1 = compute_cross1(ohm_a,j_ohm);
    SX alpha = mtimes_with_check(J_inertia_inv_sx, (cross1 + SX(tau)), "inv(J_inertia_sx)", "cross1 + tau");
    SX ohm_q = SX::zeros(4, 1);
    ohm_q(0, 0) = 0;          // Primeiro componente fixo (0)
    ohm_q(1, 0) = ohm_a(0);   // Primeiro componente de ohm_a
    ohm_q(2, 0) = ohm_a(1);   // Segundo componente de ohm_a
    ohm_q(3, 0) = ohm_a(2);   // Terceiro componente de ohm_a
    SX b = (amost * ohm_q) / 2;
    SX theta = norm_2(b); //norma de b
    SX u = b / theta;
    SX aux = u * sin(theta);

    SX cay = SX::zeros(4, 1);
    cay(0, 0) = SX(cos(theta)); // Primeiro componente: cos(theta)
    cay(1, 0) = aux(1);         // Segundo componente: primeiro valor de aux
    cay(2, 0) = aux(2);         // Terceiro componente: segundo valor de aux
    cay(3, 0) = aux(3);         // Quarto componente: terceiro valor de aux

    SX st_next_euler_e = e_a + amost * v_a;
    SX st_next_euler_v = v_a + amost * ((1/m) * mtimes_with_check(Rq, e_z, "Rq", "e_z") * torque_motor
                         - (gr * e_z));
    SX st_next_euler_ohm = ohm_a + amost*alpha;
    SX st_next_euler_q = produto_q(cay, q_a);

    SX st_next_euler = vertcat_with_check(
        {st_next_euler_e, st_next_euler_v, st_next_euler_q, st_next_euler_ohm},
        {"st_next_euler_e", "st_next_euler_v", "st_next_euler_q", "st_next_euler_ohm"}
    );


    return st_next_euler;
}

#endif
