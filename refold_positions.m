function posout=refold_positions(posin)
%Make sure that particles are inside the box
posout=posin;
posg=(posout>.5);
posl=(posout<-.5);

posout=posout-posg;
posout=posout+posl;
end
