library std;
use std.textio.all;
use std.env.all;
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use ieee.std_logic_textio.all;

use work.regionizer_data.all;


entity testbench is
--  Port ( );
end testbench;

architecture Behavioral of testbench is
    signal clk : std_logic := '0';
    signal rst : std_logic := '0';
    signal start, ready, idle, done : std_logic;
    signal newevent, newevent_out : std_logic;

    signal tk_in:  w64s(NTKSECTORS*NTKFIBERS-1 downto 0) := (others => (others => '0'));
    signal tk_out: w64s(NTKSORTED-1 downto 0) := (others => (others => '0'));
    signal calo_in:  w64s(NCALOSECTORS*NCALOFIBERS-1 downto 0) := (others => (others => '0'));
    signal calo_out: w64s(NCALOSORTED-1 downto 0) := (others => (others => '0'));
    signal mu_in:  w64s(NMUFIBERS-1 downto 0) := (others => (others => '0'));
    signal mu_out: w64s(NMUSORTED-1 downto 0) := (others => (others => '0'));

    file Fi_tk   : text open read_mode is "input-tk.txt";
    file Fi_calo : text open read_mode is "input-calo.txt";
    file Fi_mu   : text open read_mode is "input-mu.txt";
    file Fo : text open write_mode is "output-tk-vhdl_tb.txt";

begin
    clk  <= not clk after 1.25 ns;
    
    uut : entity work.full_regionizer_mux
        generic map(MU_ETA_CENTER => 460)
        port map(ap_clk => clk, 
                 ap_rst => rst, 
                 ap_start => start,
                 ap_ready => ready,
                 ap_idle =>  idle,
                 ap_done => done,
                 tracks_start => start,
                 tracks_newevent => newevent,
                 tracks_in_0_0_V => tk_in( 0),
                 tracks_in_0_1_V => tk_in( 1),
                 tracks_in_1_0_V => tk_in( 2),
                 tracks_in_1_1_V => tk_in( 3), 
                 tracks_in_2_0_V => tk_in( 4),
                 tracks_in_2_1_V => tk_in( 5),
                 tracks_in_3_0_V => tk_in( 6),
                 tracks_in_3_1_V => tk_in( 7),
                 tracks_in_4_0_V => tk_in( 8),
                 tracks_in_4_1_V => tk_in( 9), 
                 tracks_in_5_0_V => tk_in(10),
                 tracks_in_5_1_V => tk_in(11),
                 tracks_in_6_0_V => tk_in(12),
                 tracks_in_6_1_V => tk_in(13),
                 tracks_in_7_0_V => tk_in(14),
                 tracks_in_7_1_V => tk_in(15), 
                 tracks_in_8_0_V => tk_in(16),
                 tracks_in_8_1_V => tk_in(17),
                 calo_start => start,
                 calo_newevent => newevent,
                 calo_in_0_0_V => calo_in( 0),
                 calo_in_0_1_V => calo_in( 1),
                 calo_in_0_2_V => calo_in( 2),
                 calo_in_0_3_V => calo_in( 3), 
                 calo_in_1_0_V => calo_in( 4),
                 calo_in_1_1_V => calo_in( 5),
                 calo_in_1_2_V => calo_in( 6),
                 calo_in_1_3_V => calo_in( 7),
                 calo_in_2_0_V => calo_in( 8),
                 calo_in_2_1_V => calo_in( 9), 
                 calo_in_2_2_V => calo_in(10),
                 calo_in_2_3_V => calo_in(11),
                 mu_start => start,
                 mu_newevent => newevent,
                 mu_in_0_V => mu_in(0),
                 mu_in_1_V => mu_in(1),
                 tracks_out => tk_out,
                 calo_out   => calo_out,
                 mu_out     => mu_out,
                 newevent_out => newevent_out
             );
   

    runit : process 
        variable remainingEvents : integer := 5;
        variable frame : integer := 0;
        variable Li, Lo : line;
        variable itest, iobj : integer;
        variable part : particle;
        variable gpart : glbparticle;
    begin
        rst <= '1';
        wait for 5 ns;
        rst <= '0';
        start <= '0';
        tk_in <= (others => (others => '0'));
        calo_in <= (others => (others => '0'));
        mu_in <= (others => (others => '0'));
        wait until rising_edge(clk);
        while remainingEvents > 0 loop
            if not endfile(Fi_tk) then
                -- tracks
                readline(Fi_tk, Li);
                read(Li, itest);
                read(Li, iobj); if (iobj > 0) then newevent <= '1'; else newevent <= '0'; end if;
                for i in 0 to NTKSECTORS*NTKFIBERS-1  loop
                    read(Li, iobj); part.pt   := (to_signed(iobj, 16));
                    read(Li, iobj); part.eta  := (to_signed(iobj, 10));
                    read(Li, iobj); part.phi  := (to_signed(iobj, 10));
                                    part.rest := (others => '0');
                    tk_in(i) <= particle_to_w64(part);
                end loop;
                -- calo
                readline(Fi_calo, Li);
                read(Li, itest);
                read(Li, iobj); 
                for i in 0 to NCALOSECTORS*NCALOFIBERS-1  loop
                    read(Li, iobj); part.pt   := (to_signed(iobj, 16));
                    read(Li, iobj); part.eta  := (to_signed(iobj, 10));
                    read(Li, iobj); part.phi  := (to_signed(iobj, 10));
                                    part.rest := (others => '0');
                    calo_in(i) <= particle_to_w64(part);
                end loop;
                -- mu
                readline(Fi_mu, Li);
                read(Li, itest);
                read(Li, iobj); 
                for i in 0 to NMUFIBERS-1  loop
                    read(Li, iobj); gpart.pt   := (to_signed(iobj, 16));
                    read(Li, iobj); gpart.eta  := (to_signed(iobj, 12));
                    read(Li, iobj); gpart.phi  := (to_signed(iobj, 11));
                                    gpart.rest := (others => '0');
                    mu_in(i) <= glbparticle_to_w64(gpart);
                end loop;
                start <= '1';
             else
                remainingEvents := remainingEvents - 1;
                newevent <= '0';
                tk_in   <= (others => (others => '0'));
                calo_in <= (others => (others => '0'));
                mu_in   <= (others => (others => '0'));
                start <= '1';
            end if;
           -- ready to dispatch ---
            wait until rising_edge(clk);
            -- write out the output --
            write(Lo, frame, field=>5);  
            write(Lo, string'(" 1 ")); 
            write(Lo, newevent_out); 
            write(Lo, string'(" ")); 
            for i in 0 to NTKSORTED-1 loop
                part := w64_to_particle(tk_out(i));
                write(Lo, to_integer(part.pt),   field => 5); 
                write(Lo, to_integer(part.eta),  field => 5); 
                write(Lo, to_integer(part.phi),  field => 5); 
            end loop;
            for i in 0 to NCALOSORTED-1 loop
                part := w64_to_particle(calo_out(i));
                write(Lo, to_integer(part.pt),   field => 5); 
                write(Lo, to_integer(part.eta),  field => 5); 
                write(Lo, to_integer(part.phi),  field => 5); 
            end loop;
            for i in 0 to NMUSORTED-1 loop
                part := w64_to_particle(mu_out(i));
                write(Lo, to_integer(part.pt),   field => 5); 
                write(Lo, to_integer(part.eta),  field => 5); 
                write(Lo, to_integer(part.phi),  field => 5); 
            end loop;
            writeline(Fo, Lo);
            frame := frame + 1;
            --if frame >= 50 then finish(0); end if;
        end loop;
        wait for 50 ns;
        finish(0);
    end process;

    
end Behavioral;
