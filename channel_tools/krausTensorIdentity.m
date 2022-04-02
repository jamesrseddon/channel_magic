function [ new_kraus ] = krausTensorIdentity(kraus,num_qubits_A,num_qubits_B)
%KRAUSTENSORID tensor kraus operators with identity
%   Id_A \otimes kraus \otimes Id_B

num_kraus = size(kraus,3);


for kk = 1:num_kraus
    new_kraus_A = kron(eye(2^num_qubits_A),kraus(:,:,kk));
    new_kraus_B = kron(new_kraus_A,eye(2^num_qubits_B));
    new_kraus(:,:,kk) = new_kraus_B;
end
  


end

