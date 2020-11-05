library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.ipbus.all;
use work.emp_data_types.all;
use work.emp_project_decl.all;

use work.emp_device_decl.all;
use work.emp_ttc_decl.all;

use work.regionizer_data.all;

entity emp_payload is
	port(
		clk: in std_logic; -- ipbus signals
		rst: in std_logic;
		ipb_in: in ipb_wbus;
		ipb_out: out ipb_rbus;
		clk_payload: in std_logic_vector(2 downto 0);
		rst_payload: in std_logic_vector(2 downto 0);
		clk_p: in std_logic; -- data clock
		rst_loc: in std_logic_vector(N_REGION - 1 downto 0);
		clken_loc: in std_logic_vector(N_REGION - 1 downto 0);
		ctrs: in ttc_stuff_array;
		bc0: out std_logic;
		d: in ldata(4 * N_REGION - 1 downto 0); -- data in
		q: out ldata(4 * N_REGION - 1 downto 0); -- data out
		gpio: out std_logic_vector(29 downto 0); -- IO to mezzanine connector
		gpio_en: out std_logic_vector(29 downto 0) -- IO to mezzanine connector (three-state enables)
	);
		
end emp_payload;

architecture rtl of emp_payload is
        constant N_IN  : natural := NTKSECTORS*NTKFIBERS + NCALOSECTORS*NCALOFIBERS + NMUFIBERS;
        constant N_OUT : natural := NTKSORTED + NCALOSORTED + NMUSORTED;

        signal old_valid : std_logic := '0';
        signal newevent, newevent_out : std_logic;
        signal tk_in:  w64s(NTKSECTORS*NTKFIBERS-1 downto 0) := (others => (others => '0'));
        signal tk_out: w64s(NTKSORTED-1 downto 0) := (others => (others => '0'));
        signal calo_in:  w64s(NCALOSECTORS*NCALOFIBERS-1 downto 0) := (others => (others => '0'));
        signal calo_out: w64s(NCALOSORTED-1 downto 0) := (others => (others => '0'));
        signal mu_in:  w64s(NMUFIBERS-1 downto 0) := (others => (others => '0'));
        signal mu_out: w64s(NMUSORTED-1 downto 0) := (others => (others => '0'));

begin

    ipb_out <= IPB_RBUS_NULL;


    regionizer : entity work.full_regionizer_mux
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

        executor: process(clk_p)
            begin
                if rising_edge(clk_p) then
                    -- input
                    old_valid <= d(0).valid;
                    newevent  <= d(0).valid and not(old_valid);
                    for i in 0 to NTKSECTORS*NTKFIBERS-1 loop
                        tk_in(i) <= d(i).data;
                    end loop;
                    for i in 0 to NCALOSECTORS*NCALOFIBERS-1 loop
                        calo_in(i) <= d(i+NTKSECTORS*NTKFIBERS).data;
                    end loop;
                    for i in 0 to NMUFIBERS-1 loop
                        mu_in(i) <= d(i+NTKSECTORS*NTKFIBERS+NCALOSECTORS*NCALOFIBERS).data;
                    end loop;
                    -- output
                    for i in 0 to N_OUT-1 loop
                        q(i).valid  <= '1';
                        q(i).strobe <= '1';
                    end loop;
                    for i in 0 to NTKSORTED-1 loop
                        q(i).data <= tk_out(i);
                    end loop;
                    for i in 0 to NCALOSORTED-1 loop
                        q(i+NTKSORTED).data <= calo_out(i);
                    end loop;
                    for i in 0 to NMUSORTED-1 loop
                        q(i+NTKSORTED+NCALOSORTED).data <= mu_out(i);
                    end loop;
                end if;
            end process executor;

        zerofill:	
            process(clk_p) 
            begin
                if rising_edge(clk_p) then
                    for i in 4 * N_REGION - 1 downto N_OUT loop
                        q(i).data <= (others => '0');
                        q(i).valid <= '0';
                        q(i).strobe <= '1';
                    end loop;
                end if;
            end process zerofill;
    
	
	bc0 <= '0';
	
	gpio <= (others => '0');
	gpio_en <= (others => '0');

end rtl;
